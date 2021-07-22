#### RefState

using OrdinaryDiffEq
using UnPack

mutable struct ReferenceState{PS, A1}
    param_set::PS
    p0::A1
    p0_half::A1
    alpha0::A1
    alpha0_half::A1
    rho0::A1
    rho0_half::A1
    Pg::Float64
    Tg::Float64
    qtg::Float64
    sg::Float64
    function ReferenceState(Gr::Grid, param_set::PS) where {PS}

        p0 = pyzeros(Gr.nzg)
        p0_half = pyzeros(Gr.nzg)
        alpha0 = pyzeros(Gr.nzg)
        alpha0_half = pyzeros(Gr.nzg)
        rho0 = pyzeros(Gr.nzg)
        rho0_half = pyzeros(Gr.nzg)
        Pg::Float64 = 0
        Tg::Float64 = 0
        qtg::Float64 = 0
        sg::Float64 = 0
        return new{PS, typeof(p0)}(param_set, p0, p0_half, alpha0, alpha0_half, rho0, rho0_half, Pg, Tg, qtg, sg)
    end
end


"""
Initilize the reference profiles. The function is typically called from the case
specific initialization fucntion defined in Initialization.pyx
:param Gr: Grid class
:param Thermodynamics: Thermodynamics class
:param NS: StatsIO class
:param Pa:  ParallelMPI class
:return
"""
function initialize(self::ReferenceState, Gr::Grid, Stats::NetCDFIO_Stats)

    self.sg = t_to_entropy_c(self.Pg, self.Tg, self.qtg, 0.0, 0.0)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure

    function rhs(p, u, z)
        ret = eos(t_to_entropy_c, eos_first_guess_entropy, exp(p), self.qtg, self.sg)
        q_i = 0.0
        q_l = ret.ql
        T = ret.T
        return -g / (Rd * T * (1.0 - self.qtg + eps_vi * (self.qtg - q_l - q_i)))
    end

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    p0 = log(self.Pg)
    logp = p0

    p = pyzeros(Gr.nzg)
    p_half = pyzeros(Gr.nzg)

    # Perform the integration
    # TODO: replace with OrdinaryDiffEq

    z_span = (Gr.zmin, Gr.zmax)
    @show z_span
    prob = ODEProblem(rhs, logp, z_span)
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
    p_0 = [exp(sol(Gr.z[k])) for k in (Gr.gw - 1):(Gr.nzg - Gr.gw)]
    p_0_half = [exp(sol(Gr.z_half[k])) for k in (Gr.gw - 1):(Gr.nzg - Gr.gw)]

    p[(Gr.gw - 1):(Gr.nzg - Gr.gw)] .= p_0
    p_half[(Gr.gw - 1):(Gr.nzg - Gr.gw)] .= p_0_half

    # Set boundary conditions
    # TODO: re-generalize
    p[0] = p[1]
    p[end] = p[end - 1]

    p_half[0] = p_half[1]
    p_half[end] = p_half[end - 1]

    p_ = deepcopy(p)
    p_half_ = deepcopy(p_half)
    temperature = pyzeros(Gr.nzg)
    temperature_half = pyzeros(Gr.nzg)
    alpha = pyzeros(Gr.nzg)
    alpha_half = pyzeros(Gr.nzg)

    ql = pyzeros(Gr.nzg)
    qi = pyzeros(Gr.nzg)
    qv = pyzeros(Gr.nzg)

    ql_half = pyzeros(Gr.nzg)
    qi_half = pyzeros(Gr.nzg)
    qv_half = pyzeros(Gr.nzg)

    # Compute reference state thermodynamic profiles

    @inbounds for k in center_indicies(Gr)
        ret = eos(t_to_entropy_c, eos_first_guess_entropy, p_half_[k], self.qtg, self.sg)
        temperature_half[k] = ret.T
        ql_half[k] = ret.ql
        qv_half[k] = self.qtg - (ql_half[k] + qi_half[k])
        alpha_half[k] = alpha_c(p_half_[k], temperature_half[k], self.qtg, qv_half[k])
    end

    @inbounds for k in face_indicies(Gr)
        ret = eos(t_to_entropy_c, eos_first_guess_entropy, p_[k], self.qtg, self.sg)
        temperature[k] = ret.T
        ql[k] = ret.ql
        qv[k] = self.qtg - (ql[k] + qi[k])
        alpha[k] = alpha_c(p_[k], temperature[k], self.qtg, qv[k])
    end

    # Now do a sanity check to make sure that the Reference State entropy profile is uniform following
    # saturation adjustment
    local s
    @inbounds for k in xrange(Gr.nzg)
        s = t_to_entropy_c(p_half[k], temperature_half[k], self.qtg, ql_half[k], qi_half[k])
        if abs(s - self.sg) / self.sg > 0.01
            println("Error in reference profiles entropy not constant !")
            println("Likely error in saturation adjustment")
        end
    end


    self.alpha0_half = alpha_half
    self.alpha0 = alpha
    self.p0 = p_
    self.p0_half = p_half
    self.rho0 = 1.0 ./ self.alpha0
    self.rho0_half = 1.0 ./ self.alpha0_half

    # TODO: centers and faces are sliced with equal sizes,
    # they should be unequal.
    cinterior = Gr.cinterior
    finterior = Gr.finterior
    add_reference_profile(Stats, "alpha0")
    write_reference_profile(Stats, "alpha0", alpha[finterior])
    add_reference_profile(Stats, "alpha0_half")
    write_reference_profile(Stats, "alpha0_half", alpha_half[cinterior])


    add_reference_profile(Stats, "p0")
    write_reference_profile(Stats, "p0", p_[finterior])
    add_reference_profile(Stats, "p0_half")
    write_reference_profile(Stats, "p0_half", p_half[cinterior])

    add_reference_profile(Stats, "rho0")
    write_reference_profile(Stats, "rho0", 1.0 ./ alpha[finterior])
    add_reference_profile(Stats, "rho0_half")
    write_reference_profile(Stats, "rho0_half", 1.0 ./ alpha_half[cinterior])

    return
end
