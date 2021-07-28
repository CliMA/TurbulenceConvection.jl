module Cases

using Random

import ..TurbulenceConvection
using ..TurbulenceConvection: CasesBase
using ..TurbulenceConvection: set_bcs
using ..TurbulenceConvection: pow
using ..TurbulenceConvection: off_arr
using ..TurbulenceConvection: omega
using ..TurbulenceConvection: pyinterp
using ..TurbulenceConvection: eps_vi
using ..TurbulenceConvection: thetali_c
using ..TurbulenceConvection: theta_rho_c
using ..TurbulenceConvection: linspace
using ..TurbulenceConvection: eps_v
using ..TurbulenceConvection: pv_star
using ..TurbulenceConvection: fabs
using ..TurbulenceConvection: buoyancy_c
using ..TurbulenceConvection: dycoms_L
using ..TurbulenceConvection: qv_star_c
using ..TurbulenceConvection: t_to_thetali_c
using ..TurbulenceConvection: rho_c
using ..TurbulenceConvection: g
using ..TurbulenceConvection: dycoms_Rd
using ..TurbulenceConvection: dycoms_cp
using ..TurbulenceConvection: add_ts
using ..TurbulenceConvection: update
using ..TurbulenceConvection: write_ts
using ..TurbulenceConvection: initialize
using ..TurbulenceConvection: initialize_io
using ..TurbulenceConvection: io
using ..TurbulenceConvection: satadjust
using ..TurbulenceConvection: exner_c
using ..TurbulenceConvection: cpm_c
using ..TurbulenceConvection: latent_heat
using ..TurbulenceConvection: xrange
using ..TurbulenceConvection: Grid
using ..TurbulenceConvection: ReferenceState
using ..TurbulenceConvection: NetCDFIO_Stats
using ..TurbulenceConvection: GridMeanVariables
using ..TurbulenceConvection: TimeStepping
using ..TurbulenceConvection: center_indicies
using ..TurbulenceConvection: face_indicies
using ..TurbulenceConvection: real_center_indicies
using ..TurbulenceConvection: real_face_indicies

# For dispatching to inherited class
struct BaseCase end

#=
    TRMM_LBA
adopted from: "Daytime convective development over land- A model intercomparison based on LBA observations",
By Grabowski et al (2006)  Q. J. R. Meteorol. Soc. 132 317-344

    SoaresCase

Soares, P.M.M., Miranda, P.M.A., Siebesma, A.P. and Teixeira, J. (2004),
An eddy‐diffusivity/mass‐flux parametrization for dry and shallow cumulus convection.
Q.J.R. Meteorol. Soc., 130: 3365-3383. doi:10.1256/qj.03.223

    Nieuwstadt
"Nieuwstadt, F. T., Mason, P. J., Moeng, C. H., & Schumann, U. (1993).
Large-eddy simulation of the convective boundary layer: A comparison of
four computer codes. In Turbulent shear flows 8 (pp. 343-367). Springer, Berlin, Heidelberg."

    life_cycle_Tan2018
Taken from: "An extended eddy- diffusivity mass-flux scheme for unified representation of subgrid-scale turbulence and convection"
Tan, Z., Kaul, C. M., Pressel, K. G., Cohen, Y., Schneider, T., & Teixeira, J. (2018).
 Journal of Advances in Modeling Earth Systems, 10. https://doi.org/10.1002/2017MS001162

=#

struct SoaresCase end
struct BomexCase end
struct Rico end
struct TRMM_LBA end
struct ARM_SGP end
struct Nieuwstadt end
struct life_cycle_Tan2018 end

function CasesFactory(namelist, Gr, Ref)
    if namelist["meta"]["casename"] == "Soares"
        return Soares(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "Nieuwstadt"
        return Nieuwstadt(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "Bomex"
        return Bomex(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "life_cycle_Tan2018"
        return life_cycle_Tan2018(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "Rico"
        return Rico(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "TRMM_LBA"
        return TRMM_LBA(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "ARM_SGP"
        return ARM_SGP(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "GATE_III"
        return GATE_III(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "DYCOMS_RF01"
        return DYCOMS_RF01(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "GABLS"
        return GABLS(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "SP"
        return SP(namelist, Gr, Ref)
    elseif namelist["meta"]["casename"] == "DryBubble"
        return DryBubble(namelist, Gr, Ref)
    else
        error("case not recognized")
    end
    return
end

initialize_reference(self::CasesBase, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats, ::BaseCase) = nothing
initialize_profiles(self::CasesBase, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState, ::BaseCase) = nothing
initialize_surface(self::CasesBase, Gr::Grid, Ref::ReferenceState, ::BaseCase) = nothing
initialize_forcing(self::CasesBase, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables, ::BaseCase) = nothing
initialize_radiation(self::CasesBase, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables, ::BaseCase) = nothing

function TurbulenceConvection.initialize_io(self::CasesBase, Stats::NetCDFIO_Stats, ::BaseCase)
    add_ts(Stats, "Tsurface")
    add_ts(Stats, "shf")
    add_ts(Stats, "lhf")
    add_ts(Stats, "ustar")
end
function TurbulenceConvection.io(self::CasesBase, Stats::NetCDFIO_Stats, ::BaseCase)
    write_ts(Stats, "Tsurface", self.Sur.Tsurface)
    write_ts(Stats, "shf", self.Sur.shf)
    write_ts(Stats, "lhf", self.Sur.lhf)
    write_ts(Stats, "ustar", self.Sur.ustar)
end

update_surface(self::CasesBase, GMV::GridMeanVariables, TS::TimeStepping, ::BaseCase) = nothing
update_forcing(self::CasesBase, GMV::GridMeanVariables, TS::TimeStepping, ::BaseCase) = nothing
update_radiation(self::CasesBase, GMV::GridMeanVariables, TS::TimeStepping, ::BaseCase) = nothing


function Soares(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "Soares2004"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingNone}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = false
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{SoaresCase}(; casename = "SoaresCase", inversion_option, Sur, Fo, Rad)
end
function initialize_reference(self::CasesBase{SoaresCase}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1000.0 * 100.0
    Ref.qtg = 5.0e-3
    Ref.Tg = 300.0
    initialize(Ref, Gr, Stats)
end
function initialize_profiles(self::CasesBase{SoaresCase}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    theta = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0

    @inbounds for k in real_center_indicies(Gr)
        if Gr.z_half[k] <= 1350.0
            GMV.QT.values[k] = 5.0e-3 - 3.7e-4 * Gr.z_half[k] / 1000.0
            theta[k] = 300.0

        else
            GMV.QT.values[k] = 5.0e-3 - 3.7e-4 * 1.35 - 9.4e-4 * (Gr.z_half[k] - 1350.0) / 1000.0
            theta[k] = 300.0 + 2.0 * (Gr.z_half[k] - 1350.0) / 1000.0
        end
        GMV.U.values[k] = 0.01
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)

    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = theta[k]
        GMV.T.values[k] = theta[k] * exner_c(Ref.p0_half[k])
        GMV.THL.values[k] = theta[k]
    end

    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{SoaresCase}, Gr::Grid, Ref::ReferenceState)
    self.Sur.zrough = 0.16 #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    self.Sur.Tsurface = 300.0
    self.Sur.qsurface = 5.0e-3
    theta_flux = 6.0e-2
    qt_flux = 2.5e-5
    theta_surface = self.Sur.Tsurface
    self.Sur.lhf = qt_flux * Ref.rho0[Gr.gw - 1] * latent_heat(self.Sur.Tsurface) # It would be 0.0 if we follow Nieuwstadt.
    self.Sur.shf = theta_flux * cpm_c(self.Sur.qsurface) * Ref.rho0[Gr.gw - 1]
    self.Sur.ustar_fixed = false
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.bflux =
        g * (
            (theta_flux + (eps_vi - 1.0) * (theta_surface * qt_flux + self.Sur.qsurface * theta_flux)) /
            (theta_surface * (1.0 + (eps_vi - 1) * self.Sur.qsurface))
        )
    initialize(self.Sur)
end
function initialize_forcing(self::CasesBase{SoaresCase}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
end

function initialize_radiation(self::CasesBase{SoaresCase}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end



function TurbulenceConvection.initialize_io(self::CasesBase{SoaresCase}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{SoaresCase}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end
function update_surface(self::CasesBase{SoaresCase}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV)
end
function update_forcing(self::CasesBase{SoaresCase}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end

function update_radiation(self::CasesBase{SoaresCase}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

function Nieuwstadt(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "Nieuwstadt"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingNone}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = false
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{Nieuwstadt}(; casename = "Nieuwstadt", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{Nieuwstadt}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1000.0 * 100.0
    Ref.qtg = 1.0e-12 #Total water mixing ratio at TurbulenceConvection. if set to 0, alpha0, rho0, p0 are NaN (TBD)
    Ref.Tg = 300.0
    initialize(Ref, Gr, Stats)
end
function initialize_profiles(self::CasesBase{Nieuwstadt}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    theta = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0

    @inbounds for k in real_center_indicies(Gr)
        if Gr.z_half[k] <= 1350.0
            GMV.QT.values[k] = 0.0
            theta[k] = 300.0

        else
            GMV.QT.values[k] = 0.0
            theta[k] = 300.0 + 3.0 * (Gr.z_half[k] - 1350.0) / 1000.0
        end
        GMV.U.values[k] = 0.01
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)

    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = theta[k]
        GMV.T.values[k] = theta[k] * exner_c(Ref.p0_half[k])
        GMV.THL.values[k] = theta[k]
    end
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end
function initialize_surface(self::CasesBase{Nieuwstadt}, Gr::Grid, Ref::ReferenceState)
    self.Sur.zrough = 0.16 #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    self.Sur.Tsurface = 300.0
    self.Sur.qsurface = 0.0
    theta_flux = 6.0e-2
    qt_flux = 0.0
    theta_surface = self.Sur.Tsurface
    self.Sur.lhf = 0.0 # It would be 0.0 if we follow Nieuwstadt.
    self.Sur.shf = theta_flux * cpm_c(self.Sur.qsurface) * Ref.rho0[Gr.gw - 1]
    self.Sur.ustar_fixed = false
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.bflux =
        g * (
            (theta_flux + (eps_vi - 1.0) * (theta_surface * qt_flux + self.Sur.qsurface * theta_flux)) /
            (theta_surface * (1.0 + (eps_vi - 1) * self.Sur.qsurface))
        )
    initialize(self.Sur)

    return
end
function initialize_forcing(self::CasesBase{Nieuwstadt}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
end
function initialize_radiation(self::CasesBase{Nieuwstadt}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end

function TurbulenceConvection.initialize_io(self::CasesBase{Nieuwstadt}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{Nieuwstadt}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end

function update_surface(self::CasesBase{Nieuwstadt}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV)
end
function update_forcing(self::CasesBase{Nieuwstadt}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{Nieuwstadt}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

function Bomex(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "Bomex"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = true
    Fo.coriolis_param = 0.376e-4 # s^{-1}
    Fo.apply_subsidence = true
    return TurbulenceConvection.CasesBase{BomexCase}(; casename = "Bomex", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{BomexCase}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1.015e5  #Pressure at ground
    Ref.Tg = 300.4  #Temperature at ground
    Ref.qtg = 0.02245   #Total water mixing ratio at surface
    TurbulenceConvection.initialize(Ref, Gr, Stats)
end

function initialize_profiles(self::CasesBase{BomexCase}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    thetal = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0 # IC of Bomex is cloud-free

    @inbounds for k in real_center_indicies(Gr)
        #Set Thetal profile
        if Gr.z_half[k] <= 520.0
            thetal[k] = 298.7
        end
        if Gr.z_half[k] > 520.0 && Gr.z_half[k] <= 1480.0
            thetal[k] = 298.7 + (Gr.z_half[k] - 520) * (302.4 - 298.7) / (1480.0 - 520.0)
        end
        if Gr.z_half[k] > 1480.0 && Gr.z_half[k] <= 2000
            thetal[k] = 302.4 + (Gr.z_half[k] - 1480.0) * (308.2 - 302.4) / (2000.0 - 1480.0)
        end
        if Gr.z_half[k] > 2000.0
            thetal[k] = 308.2 + (Gr.z_half[k] - 2000.0) * (311.85 - 308.2) / (3000.0 - 2000.0)
        end

        #Set qt profile
        if Gr.z_half[k] <= 520
            GMV.QT.values[k] = (17.0 + (Gr.z_half[k]) * (16.3 - 17.0) / 520.0) / 1000.0
        end
        if Gr.z_half[k] > 520.0 && Gr.z_half[k] <= 1480.0
            GMV.QT.values[k] = (16.3 + (Gr.z_half[k] - 520.0) * (10.7 - 16.3) / (1480.0 - 520.0)) / 1000.0
        end
        if Gr.z_half[k] > 1480.0 && Gr.z_half[k] <= 2000.0
            GMV.QT.values[k] = (10.7 + (Gr.z_half[k] - 1480.0) * (4.2 - 10.7) / (2000.0 - 1480.0)) / 1000.0
        end
        if Gr.z_half[k] > 2000.0
            GMV.QT.values[k] = (4.2 + (Gr.z_half[k] - 2000.0) * (3.0 - 4.2) / (3000.0 - 2000.0)) / 1000.0
        end


        #Set u profile
        if Gr.z_half[k] <= 700.0
            GMV.U.values[k] = -8.75
        end
        if Gr.z_half[k] > 700.0
            GMV.U.values[k] = -8.75 + (Gr.z_half[k] - 700.0) * (-4.61 - -8.75) / (3000.0 - 700.0)
        end
    end

    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = thetal[k]
        GMV.T.values[k] = thetal[k] * exner_c(Ref.p0_half[k])
        GMV.THL.values[k] = thetal[k]
    end
    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{BomexCase}, Gr::Grid, Ref::ReferenceState)
    self.Sur.zrough = 1.0e-4 # not actually used, but initialized to reasonable value
    self.Sur.Tsurface = 299.1 * exner_c(Ref.Pg)
    self.Sur.qsurface = 22.45e-3 # kg/kg
    self.Sur.lhf = 5.2e-5 * Ref.rho0[Gr.gw - 1] * latent_heat(self.Sur.Tsurface)
    self.Sur.shf = 8.0e-3 * cpm_c(self.Sur.qsurface) * Ref.rho0[Gr.gw - 1]
    self.Sur.ustar_fixed = true
    self.Sur.ustar = 0.28 # m/s
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{BomexCase}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    @inbounds for k in real_center_indicies(Gr)
        # Geostrophic velocity profiles. vg = 0
        self.Fo.ug[k] = -10.0 + (1.8e-3) * Gr.z_half[k]
        # Set large-scale cooling
        if Gr.z_half[k] <= 1500.0
            self.Fo.dTdt[k] = (-2.0 / (3600 * 24.0)) * exner_c(Ref.p0_half[k])
        else
            self.Fo.dTdt[k] =
                (-2.0 / (3600 * 24.0) + (Gr.z_half[k] - 1500.0) * (0.0 - -2.0 / (3600 * 24.0)) / (3000.0 - 1500.0)) *
                exner_c(Ref.p0_half[k])
        end

        # Set large-scale drying
        if Gr.z_half[k] <= 300.0
            self.Fo.dqtdt[k] = -1.2e-8   #kg/(kg * s)
        end
        if Gr.z_half[k] > 300.0 && Gr.z_half[k] <= 500.0
            self.Fo.dqtdt[k] = -1.2e-8 + (Gr.z_half[k] - 300.0) * (0.0 - -1.2e-8) / (500.0 - 300.0) #kg/(kg * s)
        end

        #Set large scale subsidence
        if Gr.z_half[k] <= 1500.0
            self.Fo.subsidence[k] = 0.0 + Gr.z_half[k] * (-0.65 / 100.0 - 0.0) / (1500.0 - 0.0)
        end
        if Gr.z_half[k] > 1500.0 && Gr.z_half[k] <= 2100.0
            self.Fo.subsidence[k] = -0.65 / 100 + (Gr.z_half[k] - 1500.0) * (0.0 - -0.65 / 100.0) / (2100.0 - 1500.0)
        end
    end
    return nothing
end
function initialize_radiation(self::CasesBase{BomexCase}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end

TurbulenceConvection.initialize_io(self::CasesBase{BomexCase}, Stats::NetCDFIO_Stats) =
    initialize_io(self, Stats, BaseCase())

TurbulenceConvection.io(self::CasesBase{BomexCase}, Stats::NetCDFIO_Stats) = io(self, Stats, BaseCase())

update_surface(self::CasesBase{BomexCase}, GMV::GridMeanVariables, TS::TimeStepping) = update(self.Sur, GMV)

update_forcing(self::CasesBase{BomexCase}, GMV::GridMeanVariables, TS::TimeStepping) = update(self.Fo, GMV)
update_radiation(self::CasesBase{BomexCase}, GMV::GridMeanVariables, TS::TimeStepping) = update(self.Rad, GMV)


function life_cycle_Tan2018(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "life_cycle_Tan2018"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = true
    Fo.coriolis_param = 0.376e-4 # s^{-1}
    Fo.apply_subsidence = true
    return TurbulenceConvection.CasesBase{life_cycle_Tan2018}(;
        casename = "life_cycle_Tan2018",
        inversion_option,
        Sur,
        Fo,
        Rad,
    )
end
function initialize_reference(self::CasesBase{life_cycle_Tan2018}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1.015e5  #Pressure at ground
    Ref.Tg = 300.4  #Temperature at ground
    Ref.qtg = 0.02245   #Total water mixing ratio at surface
    initialize(Ref, Gr, Stats)
end
function initialize_profiles(self::CasesBase{life_cycle_Tan2018}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    thetal = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0 # IC of Bomex is cloud-free
    @inbounds for k in real_center_indicies(Gr)
        #Set Thetal profile
        if Gr.z_half[k] <= 520.0
            thetal[k] = 298.7
        end
        if Gr.z_half[k] > 520.0 && Gr.z_half[k] <= 1480.0
            thetal[k] = 298.7 + (Gr.z_half[k] - 520) * (302.4 - 298.7) / (1480.0 - 520.0)
        end
        if Gr.z_half[k] > 1480.0 && Gr.z_half[k] <= 2000
            thetal[k] = 302.4 + (Gr.z_half[k] - 1480.0) * (308.2 - 302.4) / (2000.0 - 1480.0)
        end
        if Gr.z_half[k] > 2000.0
            thetal[k] = 308.2 + (Gr.z_half[k] - 2000.0) * (311.85 - 308.2) / (3000.0 - 2000.0)
        end

        #Set qt profile
        if Gr.z_half[k] <= 520
            GMV.QT.values[k] = (17.0 + (Gr.z_half[k]) * (16.3 - 17.0) / 520.0) / 1000.0
        end
        if Gr.z_half[k] > 520.0 && Gr.z_half[k] <= 1480.0
            GMV.QT.values[k] = (16.3 + (Gr.z_half[k] - 520.0) * (10.7 - 16.3) / (1480.0 - 520.0)) / 1000.0
        end
        if Gr.z_half[k] > 1480.0 && Gr.z_half[k] <= 2000.0
            GMV.QT.values[k] = (10.7 + (Gr.z_half[k] - 1480.0) * (4.2 - 10.7) / (2000.0 - 1480.0)) / 1000.0
        end
        if Gr.z_half[k] > 2000.0
            GMV.QT.values[k] = (4.2 + (Gr.z_half[k] - 2000.0) * (3.0 - 4.2) / (3000.0 - 2000.0)) / 1000.0
        end


        #Set u profile
        if Gr.z_half[k] <= 700.0
            GMV.U.values[k] = -8.75
        end
        if Gr.z_half[k] > 700.0
            GMV.U.values[k] = -8.75 + (Gr.z_half[k] - 700.0) * (-4.61 - -8.75) / (3000.0 - 700.0)
        end
    end


    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = thetal[k]
        GMV.T.values[k] = thetal[k] * exner_c(Ref.p0_half[k])
        GMV.THL.values[k] = thetal[k]
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{life_cycle_Tan2018}, Gr::Grid, Ref::ReferenceState)
    self.Sur.zrough = 1.0e-4 # not actually used, but initialized to reasonable value
    self.Sur.Tsurface = 299.1 * exner_c(Ref.Pg)
    self.Sur.qsurface = 22.45e-3 # kg/kg
    self.Sur.lhf = 5.2e-5 * Ref.rho0[Gr.gw - 1] * latent_heat(self.Sur.Tsurface)
    self.Sur.shf = 8.0e-3 * cpm_c(self.Sur.qsurface) * Ref.rho0[Gr.gw - 1]
    self.lhf0 = self.Sur.lhf
    self.shf0 = self.Sur.shf
    self.Sur.ustar_fixed = true
    self.Sur.ustar = 0.28 # m/s
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.bflux = (
        g * (
            (8.0e-3 + (eps_vi - 1.0) * (299.1 * 5.2e-5 + 22.45e-3 * 8.0e-3)) /
            (299.1 * (1.0 + (eps_vi - 1) * 22.45e-3))
        )
    )
    initialize(self.Sur)
end
function initialize_forcing(self::CasesBase{life_cycle_Tan2018}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    @inbounds for k in real_center_indicies(Gr)
        # Geostrophic velocity profiles. vg = 0
        self.Fo.ug[k] = -10.0 + (1.8e-3) * Gr.z_half[k]
        # Set large-scale cooling
        if Gr.z_half[k] <= 1500.0
            self.Fo.dTdt[k] = (-2.0 / (3600 * 24.0)) * exner_c(Ref.p0_half[k])
        else
            self.Fo.dTdt[k] =
                (-2.0 / (3600 * 24.0) + (Gr.z_half[k] - 1500.0) * (0.0 - -2.0 / (3600 * 24.0)) / (3000.0 - 1500.0)) *
                exner_c(Ref.p0_half[k])
        end
        # Set large-scale drying
        if Gr.z_half[k] <= 300.0
            self.Fo.dqtdt[k] = -1.2e-8   #kg/(kg * s)
        end
        if Gr.z_half[k] > 300.0 && Gr.z_half[k] <= 500.0
            self.Fo.dqtdt[k] = -1.2e-8 + (Gr.z_half[k] - 300.0) * (0.0 - -1.2e-8) / (500.0 - 300.0) #kg/(kg * s)
        end

        #Set large scale subsidence
        if Gr.z_half[k] <= 1500.0
            self.Fo.subsidence[k] = 0.0 + Gr.z_half[k] * (-0.65 / 100.0 - 0.0) / (1500.0 - 0.0)
        end
        if Gr.z_half[k] > 1500.0 && Gr.z_half[k] <= 2100.0
            self.Fo.subsidence[k] = -0.65 / 100 + (Gr.z_half[k] - 1500.0) * (0.0 - -0.65 / 100.0) / (2100.0 - 1500.0)
        end
    end
    return nothing
end

function initialize_radiation(
    self::CasesBase{life_cycle_Tan2018},
    Gr::Grid,
    Ref::ReferenceState,
    GMV::GridMeanVariables,
)
    initialize(self.Rad, GMV)
end


function TurbulenceConvection.initialize_io(self::CasesBase{life_cycle_Tan2018}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{life_cycle_Tan2018}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end
function update_surface(self::CasesBase{life_cycle_Tan2018}, GMV::GridMeanVariables, TS::TimeStepping)
    weight = 1.0
    weight_factor = 0.01 + 0.99 * (cos(2.0 * π * TS.t / 3600.0) + 1.0) / 2.0
    weight = weight * weight_factor
    self.Sur.lhf = self.lhf0 * weight
    self.Sur.shf = self.shf0 * weight
    self.Sur.bflux = (
        g * (
            (8.0e-3 * weight + (eps_vi - 1.0) * (299.1 * 5.2e-5 * weight + 22.45e-3 * 8.0e-3 * weight)) /
            (299.1 * (1.0 + (eps_vi - 1) * 22.45e-3))
        )
    )
    update(self.Sur, GMV)
end
function update_forcing(self::CasesBase{life_cycle_Tan2018}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{life_cycle_Tan2018}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

function Rico(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "Rico"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedCoeffs; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = true
    latitude = 18.0
    Fo.coriolis_param = 2.0 * omega * sin(latitude * π / 180.0) # s^{-1}
    Fo.apply_subsidence = true
    return TurbulenceConvection.CasesBase{Rico}(; casename = "Rico", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{Rico}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1.0154e5  #Pressure at ground
    Ref.Tg = 299.8  #Temperature at ground
    pvg = pv_star(Ref.Tg)
    Ref.qtg = eps_v * pvg / (Ref.Pg - pvg)   #Total water mixing ratio at surface
    initialize(Ref, Gr, Stats)
end
function initialize_profiles(self::CasesBase{Rico}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    thetal = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0 # IC of Rico is cloud-free

    @inbounds for k in real_center_indicies(Gr)
        GMV.U.values[k] = -9.9 + 2.0e-3 * Gr.z_half[k]
        GMV.V.values[k] = -3.8
        #Set Thetal profile
        if Gr.z_half[k] <= 740.0
            thetal[k] = 297.9
        else
            thetal[k] = 297.9 + (317.0 - 297.9) / (4000.0 - 740.0) * (Gr.z_half[k] - 740.0)
        end

        #Set qt profile
        if Gr.z_half[k] <= 740.0
            GMV.QT.values[k] = (16.0 + (13.8 - 16.0) / 740.0 * Gr.z_half[k]) / 1000.0
        elseif Gr.z_half[k] > 740.0 && Gr.z_half[k] <= 3260.0
            GMV.QT.values[k] = (13.8 + (2.4 - 13.8) / (3260.0 - 740.0) * (Gr.z_half[k] - 740.0)) / 1000.0
        else
            GMV.QT.values[k] = (2.4 + (1.8 - 2.4) / (4000.0 - 3260.0) * (Gr.z_half[k] - 3260.0)) / 1000.0
        end
    end

    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = thetal[k]
        GMV.T.values[k] = thetal[k] * exner_c(Ref.p0_half[k])
        GMV.THL.values[k] = thetal[k]
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{Rico}, Gr::Grid, Ref::ReferenceState)
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.zrough = 0.00015
    self.Sur.cm = 0.001229
    self.Sur.ch = 0.001094
    self.Sur.cq = 0.001133
    # Adjust for non-IC grid spacing
    grid_adjust = (log(20.0 / self.Sur.zrough) / log(Gr.z_half[Gr.gw] / self.Sur.zrough))^2
    self.Sur.cm = self.Sur.cm * grid_adjust
    self.Sur.ch = self.Sur.ch * grid_adjust
    self.Sur.cq = self.Sur.cq * grid_adjust
    self.Sur.Tsurface = 299.8
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{Rico}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    @inbounds for k in real_center_indicies(Gr)
        # Geostrophic velocity profiles
        self.Fo.ug[k] = -9.9 + 2.0e-3 * Gr.z_half[k]
        self.Fo.vg[k] = -3.8
        # Set large-scale cooling
        self.Fo.dTdt[k] = (-2.5 / (3600.0 * 24.0)) * exner_c(Ref.p0_half[k])

        # Set large-scale moistening
        if Gr.z_half[k] <= 2980.0
            self.Fo.dqtdt[k] = (-1.0 + 1.3456 / 2980.0 * Gr.z_half[k]) / 86400.0 / 1000.0   #kg/(kg * s)
        else
            self.Fo.dqtdt[k] = 0.3456 / 86400.0 / 1000.0
        end

        #Set large scale subsidence
        if Gr.z_half[k] <= 2260.0
            self.Fo.subsidence[k] = -(0.005 / 2260.0) * Gr.z_half[k]
        else
            self.Fo.subsidence[k] = -0.005
        end
    end
    return nothing
end

function initialize_radiation(self::CasesBase{Rico}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end


function TurbulenceConvection.initialize_io(self::CasesBase{Rico}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{Rico}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end
function update_surface(self::CasesBase{Rico}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV)
end

function update_forcing(self::CasesBase{Rico}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{Rico}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

function TRMM_LBA(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "TRMM_LBA"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "thetal_maxgrad"
    Fo.apply_coriolis = false
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{TRMM_LBA}(; casename = "TRMM_LBA", inversion_option, Sur, Fo, Rad)
end
function initialize_reference(self::CasesBase{TRMM_LBA}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 991.3 * 100  #Pressure at ground
    Ref.Tg = 296.85   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    pvg = pv_star(Ref.Tg)
    Ref.qtg = eps_v * pvg / (Ref.Pg - pvg)#Total water mixing ratio at surface
    initialize(Ref, Gr, Stats)
end
function initialize_profiles(self::CasesBase{TRMM_LBA}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    p1 = TurbulenceConvection.pyzeros(Gr.nzg)

    # TRMM_LBA inputs from Grabowski et al. 2006
    #! format: off
    z_in = off_arr([0.130,  0.464,  0.573,  1.100,  1.653,  2.216,  2.760,
                     3.297,  3.824,  4.327,  4.787,  5.242,  5.686,  6.131,
                     6.578,  6.996,  7.431,  7.881,  8.300,  8.718,  9.149,
                     9.611, 10.084, 10.573, 11.008, 11.460, 11.966, 12.472,
                    12.971, 13.478, 13.971, 14.443, 14.956, 15.458, 16.019,
                    16.491, 16.961, 17.442, 17.934, 18.397, 18.851, 19.331,
                    19.809, 20.321, 20.813, 21.329, 30.000]) .* 1000 .- 130.0 #LES z is in meters

    p_in = off_arr([991.3, 954.2, 942.0, 886.9, 831.5, 778.9, 729.8,
                     684.0, 641.7, 603.2, 570.1, 538.6, 509.1, 480.4,
                     454.0, 429.6, 405.7, 382.5, 361.1, 340.9, 321.2,
                     301.2, 281.8, 263.1, 246.1, 230.1, 213.2, 197.0,
                     182.3, 167.9, 154.9, 143.0, 131.1, 119.7, 108.9,
                     100.1,  92.1,  84.6,  77.5,  71.4,  65.9,  60.7,
                      55.9,  51.3,  47.2,  43.3,  10.3]) .* 100 # LES pres is in pasc

    T_in = off_arr([23.70,  23.30,  22.57,  19.90,  16.91,  14.09,  11.13,
                      8.29,   5.38,   2.29,  -0.66,  -3.02,  -5.28,  -7.42,
                    -10.34, -12.69, -15.70, -19.21, -21.81, -24.73, -27.76,
                    -30.93, -34.62, -38.58, -42.30, -46.07, -50.03, -54.67,
                    -59.16, -63.60, -67.68, -70.77, -74.41, -77.51, -80.64,
                    -80.69, -80.00, -81.38, -81.17, -78.32, -74.77, -74.52,
                    -72.62, -70.87, -69.19, -66.90, -66.90]) .+ 273.15 # LES T is in deg K

    RH_in = off_arr([98.00,  86.00,  88.56,  87.44,  86.67,  83.67,  79.56,
                      84.78,  84.78,  89.33,  94.33,  92.00,  85.22,  77.33,
                      80.11,  66.11,  72.11,  72.67,  52.22,  54.67,  51.00,
                      43.78,  40.56,  43.11,  54.78,  46.11,  42.33,  43.22,
                      45.33,  39.78,  33.78,  28.78,  24.67,  20.67,  17.67,
                      17.11,  16.22,  14.22,  13.00,  13.00,  12.22,   9.56,
                       7.78,   5.89,   4.33,   3.00,   3.00])

    u_in = off_arr([0.00,   0.81,   1.17,   3.44,   3.53,   3.88,   4.09,
                     3.97,   1.22,   0.16,  -1.22,  -1.72,  -2.77,  -2.65,
                    -0.64,  -0.07,  -1.90,  -2.70,  -2.99,  -3.66,  -5.05,
                    -6.64,  -4.74,  -5.30,  -6.07,  -4.26,  -7.52,  -8.88,
                    -9.00,  -7.77,  -5.37,  -3.88,  -1.15,  -2.36,  -9.20,
                    -8.01,  -5.68,  -8.83, -14.51, -15.55, -15.36, -17.67,
                   -17.82, -18.94, -15.92, -15.32, -15.32])

    v_in = off_arr([-0.40,  -3.51,  -3.88,  -4.77,  -5.28,  -5.85,  -5.60,
                     -2.67,  -1.47,   0.57,   0.89,  -0.08,   1.11,   2.15,
                      3.12,   3.22,   3.34,   1.91,   1.15,   1.01,  -0.57,
                     -0.67,   0.31,   2.97,   2.32,   2.66,   4.79,   3.40,
                      3.14,   3.93,   7.57,   2.58,   2.50,   6.44,   6.84,
                      0.19,  -2.20,  -3.60,   0.56,   6.68,   9.41,   7.03,
                      5.32,   1.14,  -0.65,   5.27,   5.27])
    #! format: on
    # interpolate to the model grid-points

    p1 = pyinterp(Gr.z_half, z_in, p_in)
    GMV.U.values .= pyinterp(Gr.z_half, z_in, u_in)
    GMV.V.values .= pyinterp(Gr.z_half, z_in, v_in)

    # get the entropy from RH, p, T
    RH = TurbulenceConvection.pyzeros(Gr.nzg)
    z_half_in = off_arr(Gr.z_half[(Gr.gw):(Gr.nzg - Gr.gw)])
    RH[(Gr.gw):(Gr.nzg - Gr.gw)] = pyinterp(z_half_in, z_in, RH_in)
    RH[0] = RH[3]
    RH[1] = RH[2]
    RH[Gr.nzg - Gr.gw + 1] = RH[Gr.nzg - Gr.gw - 1]

    T = TurbulenceConvection.pyzeros(Gr.nzg)
    T[(Gr.gw):(Gr.nzg - Gr.gw)] = pyinterp(z_half_in, z_in, T_in)
    GMV.T.values .= T
    theta_rho = RH * 0.0
    epsi = 287.1 / 461.5

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.T, Gr)


    @inbounds for k in real_center_indicies(Gr)
        PV_star = pv_star(GMV.T.values[k])
        qv_star = PV_star * epsi / (p1[k] - PV_star + epsi * PV_star * RH[k] / 100.0) # eq. 37 in pressel et al and the def of RH
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        GMV.QT.values[k] = qv_star * RH[k] / 100.0
        GMV.H.values[k] =
            thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], 0.0, 0.0, latent_heat(GMV.T.values[k]))

        GMV.THL.values[k] =
            thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], 0.0, 0.0, latent_heat(GMV.T.values[k]))
        theta_rho[k] = theta_rho_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], qv)
    end

    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{TRMM_LBA}, Gr::Grid, Ref::ReferenceState)
    #self.Sur.zrough = 1.0e-4 # not actually used, but initialized to reasonable value
    self.Sur.Tsurface = (273.15 + 23) * exner_c(Ref.Pg)
    self.Sur.qsurface = 22.45e-3 # kg/kg
    self.Sur.lhf = 5.2e-5 * Ref.rho0[Gr.gw - 1] * latent_heat(self.Sur.Tsurface)
    self.Sur.shf = 8.0e-3 * cpm_c(self.Sur.qsurface) * Ref.rho0[Gr.gw - 1]
    self.Sur.ustar_fixed = true
    self.Sur.ustar = 0.28 # this is taken from Bomex -- better option is to approximate from LES tke above the surface
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{TRMM_LBA}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    self.Fo.dTdt = TurbulenceConvection.pyzeros(Gr.nzg)
    self.rad_time = linspace(10, 360; num = 36) .* 60
    #! format: off
    z_in         = off_arr([42.5, 200.92, 456.28, 743, 1061.08, 1410.52, 1791.32, 2203.48, 2647,3121.88, 3628.12,
                             4165.72, 4734.68, 5335, 5966.68, 6629.72, 7324.12,
                             8049.88, 8807, 9595.48, 10415.32, 11266.52, 12149.08, 13063, 14008.28,
                             14984.92, 15992.92, 17032.28, 18103, 19205.08, 20338.52, 21503.32, 22699.48])
    rad_in   = [[-1.386, -1.927, -2.089, -1.969, -1.805, -1.585, -1.406, -1.317, -1.188, -1.106, -1.103, -1.025,
                          -0.955, -1.045, -1.144, -1.119, -1.068, -1.092, -1.196, -1.253, -1.266, -1.306,  -0.95,  0.122,
                           0.255,  0.258,  0.322,  0.135,      0,      0,      0,      0,      0],
                         [ -1.23, -1.824, -2.011, -1.895, -1.729, -1.508, -1.331, -1.241, -1.109, -1.024, -1.018,  -0.94,
                          -0.867, -0.953, -1.046, -1.018, -0.972, -1.006, -1.119, -1.187, -1.209, -1.259, -0.919,  0.122,
                           0.264,  0.262,  0.326,  0.137,      0,      0,      0,      0,     0],
                         [-1.043, -1.692, -1.906, -1.796,  -1.63,  -1.41, -1.233, -1.142,  -1.01,  -0.92, -0.911, -0.829,
                          -0.754, -0.837, -0.923,  -0.89, -0.847, -0.895, -1.021, -1.101, -1.138, -1.201,  -0.88,  0.131,
                           0.286,  0.259,  0.332,   0.14,      0,      0,      0,      0,      0],
                         [-0.944, -1.613, -1.832,  -1.72, -1.555, -1.339, -1.163, -1.068, -0.935, -0.846, -0.835,  -0.75,
                          -0.673, -0.751, -0.833, -0.798,  -0.76, -0.817, -0.952, -1.042, -1.088, -1.159, -0.853,  0.138,
                           0.291,  0.265,  0.348,  0.136,      0,      0,      0,      0,      0],
                         [-0.833, -1.526, -1.757, -1.648, -1.485,  -1.27, -1.093, -0.998, -0.867, -0.778, -0.761, -0.672,
                          -0.594, -0.671, -0.748, -0.709, -0.676, -0.742, -0.887, -0.986, -1.041, -1.119, -0.825,  0.143,
                           0.296,  0.271,  0.351,  0.138,      0,      0,      0,      0,      0],
                         [-0.719, -1.425, -1.657,  -1.55, -1.392, -1.179, -1.003, -0.909, -0.778, -0.688, -0.667, -0.573,
                          -0.492, -0.566, -0.639, -0.596, -0.568, -0.647, -0.804, -0.914, -0.981,  -1.07, -0.793,  0.151,
                           0.303,  0.279,  0.355,  0.141,      0,      0,      0,      0,      0],
                         [-0.724, -1.374, -1.585, -1.482, -1.328, -1.116, -0.936, -0.842, -0.715, -0.624, -0.598, -0.503,
                          -0.421, -0.494, -0.561, -0.514,  -0.49,  -0.58, -0.745, -0.863, -0.938, -1.035, -0.764,  0.171,
                           0.291,  0.284,  0.358,  0.144,      0,      0,      0,      0,      0],
                         [-0.587,  -1.28, -1.513, -1.416, -1.264, -1.052, -0.874, -0.781, -0.655, -0.561, -0.532, -0.436,
                          -0.354, -0.424, -0.485, -0.435, -0.417, -0.517, -0.691, -0.817, -0.898,     -1,  -0.74,  0.176,
                           0.297,  0.289,   0.36,  0.146,      0,      0,      0,      0,      0],
                         [-0.506, -1.194, -1.426, -1.332, -1.182, -0.972, -0.795, -0.704, -0.578,  -0.48, -0.445, -0.347,
                          -0.267, -0.336, -0.391, -0.337, -0.325, -0.436,  -0.62, -0.756, -0.847,  -0.96, -0.714,   0.18,
                           0.305,  0.317,  0.348,  0.158,      0,      0,      0,      0,      0],
                         [-0.472,  -1.14, -1.364, -1.271, -1.123, -0.914, -0.738, -0.649, -0.522, -0.422, -0.386, -0.287,
                          -0.207, -0.273, -0.322, -0.267,  -0.26, -0.379, -0.569, -0.712, -0.811, -0.931, -0.696,  0.183,
                           0.311,   0.32,  0.351,   0.16,      0,      0,      0,      0,     0],
                         [-0.448, -1.091, -1.305, -1.214, -1.068, -0.858, -0.682, -0.594, -0.469, -0.368, -0.329, -0.229,
                          -0.149, -0.213, -0.257,   -0.2, -0.199, -0.327, -0.523, -0.668, -0.774, -0.903, -0.678,  0.186,
                           0.315,  0.323,  0.355,  0.162,      0,      0,      0,      0,      0],
                         [-0.405, -1.025, -1.228, -1.139, -0.996, -0.789, -0.615, -0.527, -0.402,   -0.3, -0.256, -0.156,
                          -0.077, -0.136, -0.173, -0.115, -0.121, -0.259, -0.463, -0.617, -0.732, -0.869, -0.656,   0.19,
                           0.322,  0.326,  0.359,  0.164,      0,      0,      0,      0,      0],
                         [-0.391, -0.983, -1.174, -1.085, -0.945, -0.739, -0.566, -0.478, -0.354, -0.251, -0.205, -0.105,
                          -0.027, -0.082, -0.114, -0.056, -0.069, -0.213,  -0.42, -0.579, -0.699,  -0.84, -0.642,  0.173,
                           0.327,  0.329,  0.362,  0.165,      0,      0,      0,      0,      0],
                         [-0.385, -0.946, -1.121, -1.032, -0.898, -0.695, -0.523, -0.434, -0.307, -0.203, -0.157, -0.057,
                           0.021, -0.031, -0.059, -0.001, -0.018, -0.168, -0.381, -0.546, -0.672, -0.819, -0.629,  0.176,
                           0.332,  0.332,  0.364,  0.166,      0,      0,      0,      0,      0],
                         [-0.383, -0.904, -1.063, -0.972, -0.834, -0.632, -0.464, -0.378, -0.251, -0.144, -0.096,  0.001,
                           0.079,  0.032,  0.011,  0.069,  0.044, -0.113, -0.332, -0.504, -0.637, -0.791, -0.611,  0.181,
                           0.338,  0.335,  0.367,  0.167,      0,      0,      0,      0,      0],
                         [-0.391, -0.873, -1.016, -0.929, -0.794, -0.591, -0.423, -0.337, -0.212, -0.104, -0.056,  0.043,
                           0.121,  0.077,  0.058,  0.117,  0.088, -0.075, -0.298, -0.475, -0.613, -0.772, -0.599,  0.183,
                           0.342,  0.337,   0.37,  0.168,      0,      0,      0,      0,      0],
                         [-0.359, -0.836, -0.976, -0.888, -0.755, -0.554, -0.386,   -0.3, -0.175, -0.067, -0.018,  0.081,
                            0.16,  0.119,  0.103,  0.161,  0.129, -0.039, -0.266, -0.448, -0.591, -0.755, -0.587,  0.187,
                           0.345,  0.339,  0.372,  0.169,      0,      0,      0,      0,     0],
                         [-0.328, -0.792, -0.928, -0.842, -0.709, -0.508, -0.341, -0.256, -0.131, -0.022,  0.029,  0.128,
                           0.208,   0.17,  0.158,  0.216,  0.179,  0.005, -0.228, -0.415, -0.564, -0.733, -0.573,   0.19,
                           0.384,  0.313,  0.375,   0.17,      0,      0,      0,      0,      0],
                         [-0.324, -0.767, -0.893, -0.807, -0.676, -0.476,  -0.31, -0.225, -0.101,  0.008,   0.06,  0.159,
                           0.239,  0.204,  0.195,  0.252,  0.212,  0.034, -0.203, -0.394, -0.546, -0.719, -0.564,  0.192,
                           0.386,  0.315,  0.377,  0.171,      0,      0,      0,      0,      0],
                         [ -0.31,  -0.74,  -0.86, -0.775, -0.647, -0.449, -0.283, -0.197, -0.073,  0.036,  0.089,  0.188,
                           0.269,  0.235,  0.229,  0.285,  0.242,  0.061, -0.179, -0.374,  -0.53, -0.706, -0.556,  0.194,
                           0.388,  0.317,  0.402,  0.158,      0,      0,      0,      0,      0],
                         [-0.244, -0.694, -0.818,  -0.73, -0.605, -0.415, -0.252, -0.163, -0.037,  0.072,  0.122,   0.22,
                           0.303,  0.273,  0.269,  0.324,  0.277,  0.093, -0.152,  -0.35,  -0.51, -0.691, -0.546,  0.196,
                           0.39,   0.32,  0.403,  0.159,      0,      0,      0,      0,      0],
                         [-0.284, -0.701, -0.803, -0.701, -0.568, -0.381, -0.225, -0.142, -0.017,  0.092,  0.143,  0.242,
                           0.325,  0.298,  0.295,   0.35,    0.3,  0.112, -0.134, -0.334, -0.497,  -0.68,  -0.54,  0.198,
                           0.392,  0.321,  0.404,   0.16,      0,      0,      0,      0,      0],
                         [-0.281, -0.686, -0.783,  -0.68, -0.547, -0.359, -0.202, -0.119,  0.005,  0.112,  0.163,  0.261,
                           0.345,  0.321,  0.319,  0.371,  0.319,   0.13, -0.118, -0.321, -0.486, -0.671, -0.534,  0.199,
                           0.393,  0.323,  0.405,  0.161,      0,      0,      0,      0,      0],
                         [-0.269, -0.667,  -0.76, -0.655, -0.522, -0.336, -0.181, -0.096,  0.029,  0.136,  0.188,  0.286,
                            0.37,  0.346,  0.345,  0.396,  0.342,   0.15, -0.102, -0.307, -0.473, -0.661, -0.528,    0.2,
                           0.393,  0.324,  0.405,  0.162,      0,      0,      0,      0,      0],
                         [-0.255, -0.653, -0.747, -0.643, -0.511, -0.325, -0.169, -0.082,  0.042,  0.149,  0.204,  0.304,
                           0.388,  0.363,  0.36 ,  0.409,  0.354,  0.164, -0.085, -0.289, -0.457, -0.649, -0.523,  0.193,
                           0.394,  0.326,  0.406,  0.162,      0,      0,      0,      0,      0],
                         [-0.265,  -0.65, -0.739, -0.634,   -0.5, -0.314, -0.159, -0.072,  0.052,  0.159,  0.215,  0.316,
                           0.398,  0.374,  0.374,  0.424,   0.37,  0.181, -0.065, -0.265, -0.429, -0.627, -0.519,   0.18,
                           0.394,  0.326,  0.406,  0.162,      0,      0,      0,      0,      0],
                         [-0.276, -0.647, -0.731, -0.626, -0.492, -0.307, -0.152, -0.064,  0.058,  0.166,  0.227,  0.329,
                           0.411,  0.389,   0.39,  0.441,  0.389,  0.207, -0.032, -0.228, -0.394, -0.596, -0.494,  0.194,
                           0.376,  0.326,  0.406,  0.162,      0,      0,      0,      0,      0],
                         [-0.271, -0.646,  -0.73, -0.625, -0.489, -0.303, -0.149, -0.061,  0.062,  0.169,  0.229,  0.332,
                           0.412,  0.388,  0.389,  0.439,  0.387,  0.206, -0.028, -0.209, -0.347, -0.524, -0.435,  0.195,
                           0.381,  0.313,  0.405,  0.162,      0,      0,      0,      0,      0],
                         [-0.267, -0.647, -0.734, -0.628,  -0.49, -0.304, -0.151, -0.062,  0.061,  0.168,  0.229,  0.329,
                           0.408,  0.385,  0.388,  0.438,  0.386,  0.206, -0.024, -0.194, -0.319,  -0.48,  -0.36,  0.318,
                           0.405,  0.335,  0.394,  0.162,      0,      0,      0,      0,      0],
                         [-0.274, -0.656, -0.745,  -0.64,   -0.5, -0.313, -0.158, -0.068,  0.054,  0.161,  0.223,  0.325,
                           0.402,  0.379,  0.384,  0.438,  0.392,  0.221,  0.001, -0.164, -0.278, -0.415, -0.264,  0.445,
                           0.402,  0.304,  0.389,  0.157,      0,      0,      0,      0,      0],
                         [-0.289, -0.666, -0.753, -0.648, -0.508,  -0.32, -0.164, -0.073,  0.049,  0.156,   0.22,  0.321,
                           0.397,  0.374,  0.377,   0.43,  0.387,  0.224,  0.014, -0.139, -0.236, -0.359, -0.211,  0.475,
                             0.4,  0.308,  0.375,  0.155,      0,      0,      0,      0,      0],
                         [-0.302, -0.678, -0.765, -0.659, -0.517, -0.329, -0.176, -0.085,  0.038,  0.145,  0.208,   0.31,
                           0.386,  0.362,  0.366,  0.421,  0.381,  0.224,  0.022, -0.119, -0.201,   -0.3, -0.129,  0.572,
                           0.419,  0.265,  0.364,  0.154,      0,      0,      0,      0,      0],
                         [-0.314, -0.696, -0.786, -0.681, -0.539, -0.349, -0.196, -0.105,  0.019,  0.127,  0.189,  0.289,
                           0.364,   0.34,  0.346,  0.403,   0.37,  0.222,  0.036, -0.081, -0.133, -0.205, -0.021,  0.674,
                           0.383,  0.237,  0.359,  0.151,      0,      0,      0,      0,      0],
                         [-0.341, -0.719, -0.807, -0.702, -0.558, -0.367, -0.211,  -0.12,  0.003,  0.111,  0.175,  0.277,
                           0.351,  0.325,  0.331,   0.39,   0.36,  0.221,  0.048, -0.046, -0.074, -0.139,  0.038,  0.726,
                           0.429,  0.215,  0.347,  0.151,      0,      0,      0,      0,      0],
                         [ -0.35, -0.737, -0.829, -0.724, -0.577, -0.385, -0.229, -0.136, -0.011,  0.098,  0.163,  0.266,
                           0.338,   0.31,  0.316,  0.378,  0.354,  0.221,  0.062, -0.009, -0.012, -0.063,  0.119,  0.811,
                           0.319,  0.201,  0.343,  0.148,      0,      0,      0,      0,      0],
                         [-0.344,  -0.75, -0.856, -0.757, -0.607, -0.409,  -0.25, -0.156, -0.033,  0.076,  0.143,  0.246,
                           0.316,  0.287,  0.293,  0.361,  0.345,  0.225,  0.082,  0.035,  0.071,  0.046,  0.172,  0.708,
                           0.255,   0.21,  0.325,  0.146,      0,      0,      0,      0,      0]] ./ 86400

    #! format: on
    # TODO: check translation
    rad_in = reduce(vcat, rad_in')
    rad_in = off_arr(rad_in)
    A = hcat(map(xrange(0, 36)) do tt
        a = Gr.z_half
        b = z_in
        c = reshape(rad_in[tt, :], size(rad_in, 2))
        c = off_arr(c)
        pyinterp(a, b, c)
    end...)
    A = off_arr(A')

    self.rad = A # store matrix in self
    ind1 = Int(trunc(10.0 / 600.0))
    ind2 = Int(ceil(10.0 / 600.0)) - 1
    @inbounds for k in face_indicies(Gr)
        if 10 % 600.0 == 0
            self.Fo.dTdt[k] = self.rad[ind1, k]
        else
            self.Fo.dTdt[k] =
                (self.rad[ind2, k] - self.rad[ind1, k]) / (self.rad_time[ind2 + 1] - self.rad_time[ind1 + 1]) * (10.0) +
                self.rad[ind1, k]
        end
    end
    return nothing
end
function initialize_radiation(self::CasesBase{TRMM_LBA}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end


function TurbulenceConvection.initialize_io(self::CasesBase{TRMM_LBA}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{TRMM_LBA}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end

function update_surface(self::CasesBase{TRMM_LBA}, GMV::GridMeanVariables, TS::TimeStepping)
    self.Sur.lhf = 554.0 * pow(max(0, cos(π / 2 * ((5.25 * 3600.0 - TS.t) / 5.25 / 3600.0))), 1.3)
    self.Sur.shf = 270.0 * pow(max(0, cos(π / 2 * ((5.25 * 3600.0 - TS.t) / 5.25 / 3600.0))), 1.5)
    update(self.Sur, GMV)
    # fix momentum fluxes to zero as they are not used in the paper
    self.Sur.rho_uflux = 0.0
    self.Sur.rho_vflux = 0.0
end

function update_forcing(self::CasesBase{TRMM_LBA}, GMV::GridMeanVariables, TS::TimeStepping)

    Gr = self.Fo.Gr
    ind2 = Int(ceil(TS.t / 600.0))
    ind1 = Int(trunc(TS.t / 600.0))
    @inbounds for k in face_indicies(Gr)
        if Gr.z_half[k] >= 22699.48
            self.Fo.dTdt[k] = 0.0
        else
            if TS.t < 600.0 # first 10 min use the radiative forcing of t=10min (as in the paper)
                self.Fo.dTdt[k] = self.rad[0, k]
            elseif TS.t < 21600.0 && ind2 < 36
                if TS.t % 600.0 == 0
                    self.Fo.dTdt[k] = self.rad[ind1, k]
                else
                    self.Fo.dTdt[k] =
                        (self.rad[ind2, k] - self.rad[ind1, k]) / (self.rad_time[ind2] - self.rad_time[ind1]) *
                        (TS.t - self.rad_time[ind1]) + self.rad[ind1, k]
                end
            else
                self.Fo.dTdt[k] = self.rad[35, k]
            end
        end
    end

    update(self.Fo, GMV)

end
function update_radiation(self::CasesBase{TRMM_LBA}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

# adopted from: "Large-eddy simulation of the diurnal cycle of shallow cumulus convection over land",
# By Brown et al. (2002)  Q. J. R. Meteorol. Soc. 128, 1075-1093
function ARM_SGP(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "ARM_SGP"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "thetal_maxgrad"
    Fo.apply_coriolis = true
    Fo.coriolis_param = 8.5e-5
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{ARM_SGP}(; casename = "ARM_SGP", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{ARM_SGP}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 970.0 * 100 #Pressure at ground
    Ref.Tg = 299.0   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    Ref.qtg = 15.2 / 1000#Total water mixing ratio at surface
    initialize(Ref, Gr, Stats)
end

function initialize_profiles(self::CasesBase{ARM_SGP}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    p1 = TurbulenceConvection.pyzeros(Gr.nzg)

    # ARM_SGP inputs
    #! format: off
    z_in = off_arr([0.0, 50.0, 350.0, 650.0, 700.0, 1300.0, 2500.0, 5500.0 ]) #LES z is in meters
    Theta_in = off_arr([299.0, 301.5, 302.5, 303.53, 303.7, 307.13, 314.0, 343.2]) # K
    r_in = off_arr([15.2,15.17,14.98,14.8,14.7,13.5,3.0,3.0])/1000 # qt should be in kg/kg
    #! format: on
    qt_in = r_in ./ (1 .+ r_in)

    # interpolate to the model grid-points
    Theta = pyinterp(Gr.z_half, z_in, Theta_in)
    qt = pyinterp(Gr.z_half, z_in, qt_in)


    @inbounds for k in real_center_indicies(Gr)
        GMV.U.values[k] = 10.0
        GMV.QT.values[k] = qt[k]
        GMV.T.values[k] = Theta[k] * exner_c(Ref.p0_half[k])
        GMV.H.values[k] =
            thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], 0.0, 0.0, latent_heat(GMV.T.values[k]))

        GMV.THL.values[k] =
            thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], 0.0, 0.0, latent_heat(GMV.T.values[k]))
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{ARM_SGP}, Gr::Grid, Ref::ReferenceState)
    self.Sur.Tsurface = 299.0 * exner_c(Ref.Pg)
    self.Sur.qsurface = 15.2e-3 # kg/kg
    self.Sur.lhf = 5.0
    self.Sur.shf = -30.0
    self.Sur.ustar_fixed = true
    self.Sur.ustar = 0.28 # this is taken from Bomex -- better option is to approximate from LES tke above the surface
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{ARM_SGP}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    @inbounds for k in center_indicies(Gr)
        self.Fo.ug[k] = 10.0
        self.Fo.vg[k] = 0.0
    end
    return nothing
end

function initialize_radiation(self::CasesBase{ARM_SGP}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end

function TurbulenceConvection.initialize_io(self::CasesBase{ARM_SGP}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{ARM_SGP}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end

function update_surface(self::CasesBase{ARM_SGP}, GMV::GridMeanVariables, TS::TimeStepping)
    t_Sur_in = off_arr([0.0, 4.0, 6.5, 7.5, 10.0, 12.5, 14.5]) .* 3600 #LES time is in sec
    SH = off_arr([-30.0, 90.0, 140.0, 140.0, 100.0, -10, -10]) # W/m^2
    LH = off_arr([5.0, 250.0, 450.0, 500.0, 420.0, 180.0, 0.0]) # W/m^2
    self.Sur.shf = pyinterp(off_arr([TS.t]), t_Sur_in, SH)[0]
    self.Sur.lhf = pyinterp(off_arr([TS.t]), t_Sur_in, LH)[0]
    # if fluxes vanish bflux vanish and wstar and obukov length are NaNs
    ## CK +++ I commented out the lines below as I don"t think this is how we want to fix things!
    # if self.Sur.shf < 1.0
    #     self.Sur.shf = 1.0
    # if self.Sur.lhf < 1.0
    #     self.Sur.lhf = 1.0
    #+++++++++
    update(self.Sur, GMV)
    # fix momentum fluxes to zero as they are not used in the paper
    self.Sur.rho_uflux = 0.0
    self.Sur.rho_vflux = 0.0
end

function update_forcing(self::CasesBase{ARM_SGP}, GMV::GridMeanVariables, TS::TimeStepping)
    t_in = off_arr([0.0, 3.0, 6.0, 9.0, 12.0, 14.5]) .* 3600.0 #LES time is in sec
    AT_in = off_arr([0.0, 0.0, 0.0, -0.08, -0.016, -0.016]) ./ 3600.0 # Advective forcing for theta [K/h] converted to [K/sec]
    RT_in = off_arr([-0.125, 0.0, 0.0, 0.0, 0.0, -0.1]) ./ 3600.0  # Radiative forcing for theta [K/h] converted to [K/sec]
    Rqt_in = off_arr([0.08, 0.02, 0.04, -0.1, -0.16, -0.3]) ./ 1000.0 ./ 3600.0 # Radiative forcing for qt converted to [kg/kg/sec]
    dTdt = pyinterp(off_arr([TS.t]), t_in, AT_in)[0] + pyinterp(off_arr([TS.t]), t_in, RT_in)[0]
    dqtdt = pyinterp(off_arr([TS.t]), t_in, Rqt_in)[0]
    Gr = self.Fo.Gr
    @inbounds for k in center_indicies(Gr)
        if Gr.z_half[k] <= 1000.0
            self.Fo.dTdt[k] = dTdt
            self.Fo.dqtdt[k] = dqtdt * exner_c(self.Fo.Ref.p0_half[k])
        elseif Gr.z_half[k] > 1000.0 && Gr.z_half[k] <= 2000.0
            self.Fo.dTdt[k] = dTdt * (1 - (Gr.z_half[k] - 1000.0) / 1000.0)
            self.Fo.dqtdt[k] = dqtdt * exner_c(self.Fo.Ref.p0_half[k]) * (1 - (Gr.z_half[k] - 1000.0) / 1000.0)
        end
    end
    update(self.Fo, GMV)
end

function update_radiation(self::CasesBase{ARM_SGP}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

# adopted from: "Large eddy simulation of Maritime Deep Tropical Convection",
# By Khairoutdinov et al (2009)  JAMES, vol. 1, article #15
function GATE_III(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "GATE_III"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedCoeffs; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "thetal_maxgrad"
    Fo.apply_subsidence = false
    Fo.apply_coriolis = false
    return TurbulenceConvection.CasesBase{GATE_III}(; casename = "GATE_III", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{GATE_III}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1013.0 * 100  #Pressure at ground
    Ref.Tg = 299.184   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    Ref.qtg = 16.5 / 1000#Total water mixing ratio at surface
    initialize(Ref, Gr, Stats)
end

function initialize_profiles(self::CasesBase{GATE_III}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    qt = TurbulenceConvection.pyzeros(Gr.nzg)
    T = TurbulenceConvection.pyzeros(Gr.nzg) # Gr.nzg = Gr.nz + 2*Gr.gw
    U = TurbulenceConvection.pyzeros(Gr.nzg)
    theta_rho = TurbulenceConvection.pyzeros(Gr.nzg)

    # GATE_III inputs - I extended them to z=22 km
    #! format: off
    z_in  = off_arr([ 0.0,   0.5,  1.0,  1.5,  2.0,   2.5,    3.0,   3.5,   4.0,   4.5,   5.0,  5.5,  6.0,  6.5,
                       7.0, 7.5, 8.0,  8.5,   9.0,   9.5, 10.0,   10.5,   11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0,
                       14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 27.0]) * 1000.0 #z is in meters
    r_in = off_arr([16.5,  16.5, 13.5, 12.0, 10.0,   8.7,    7.1,   6.1,   5.2,   4.5,   3.6,  3.0,  2.3, 1.75, 1.3,
                     0.9, 0.5, 0.25, 0.125, 0.065, 0.003, 0.0015, 0.0007,  0.0003,  0.0001,  0.0001,  0.0001,  0.0001,
                     0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001, 0.0001])/1000 # mixing ratio should be in kg/kg
    U_in  = off_arr([  -1, -1.75, -2.5, -3.6, -6.0, -8.75, -11.75, -13.0, -13.1, -12.1, -11.0, -8.5, -5.0, -2.6, 0.0,
                        0.5, 0.4,  0.3,   0.0,  -1.0, -2.5,   -3.5,   -4.5, -4.8, -5.0, -3.5, -2.0, -1.0, -1.0, -1.0,
                        -1.5, -2.0, -2.5, -2.6, -2.7, -3.0, -3.0, -3.0])# [m/s]
    #! format: on
    qt_in = r_in ./ (1 .+ r_in) # convert mixing ratio to specific humidity

    # temperature is taken from a different input plot at different z levels
    #! format: off
    T_in = off_arr([299.184, 294.836, 294.261, 288.773, 276.698, 265.004, 253.930, 243.662, 227.674, 214.266, 207.757, 201.973, 198.278, 197.414, 198.110, 198.110])
    z_T_in = off_arr([0.0, 0.492, 0.700, 1.698, 3.928, 6.039, 7.795, 9.137, 11.055, 12.645, 13.521, 14.486, 15.448, 16.436, 17.293, 22.0])*1000.0 # for km
    #! format: on

    # interpolate to the model grid-points
    T = pyinterp(Gr.z_half, z_T_in, T_in) # interpolate to ref pressure level
    qt = pyinterp(Gr.z_half, z_in, qt_in)
    U = pyinterp(Gr.z_half, z_in, U_in)


    @inbounds for k in real_center_indicies(Gr)
        GMV.QT.values[k] = qt[k]
        GMV.T.values[k] = T[k]
        GMV.U.values[k] = U[k]

        GMV.H.values[k] =
            thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], 0.0, 0.0, latent_heat(GMV.T.values[k]))

        GMV.THL.values[k] =
            thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], 0.0, 0.0, latent_heat(GMV.T.values[k]))
    end
    set_bcs(GMV.U, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.T, Gr)
    set_bcs(GMV.H, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{GATE_III}, Gr::Grid, Ref::ReferenceState)
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.qsurface = 16.5 / 1000.0 # kg/kg
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.cm = 0.0012
    self.Sur.ch = 0.0034337
    self.Sur.cq = 0.0034337
    self.Sur.Tsurface = 299.184
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{GATE_III}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    #LES z is in meters
    #! format: off
    z_in     = off_arr([ 0.0,   0.5,  1.0,  1.5,   2.0,   2.5,    3.0,   3.5,   4.0,   4.5,   5.0,   5.5,   6.0,
                          6.5,  7.0,  7.5,   8.0,  8.5,   9.0,  9.5,  10.0,  10.5,  11.0,    11.5,   12.0, 12.5,
                          13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0]) * 1000.0
    u_in     = off_arr([  -1, -1.75, -2.5, -3.6,  -6.0, -8.75, -11.75, -12.9, -13.1, -12.1, -11.0,  -8.5,  -5.0,
                           -2.6,  0.0,  0.5,   0.4,  0.3,   0.0, -1.0,  -3.0,  -3.5,  -4.5,    -4.6,   -5.0, -3.5,
                           -2.0, -1.0, -1.0, -1.0, -1.5, -2.0, -2.5, -2.6, -2.7, -3.0, -3.0])
    # Radiative forcing for T [K/d] converted to [K/sec]
    RAD_in   = off_arr([-2.9,  -1.1, -0.8, -1.1, -1.25, -1.35,   -1.4,  -1.4, -1.44, -1.52,  -1.6, -1.54, -1.49,
                         -1.43, -1.36, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05,  -1.0,  -0.95,   -0.9,  -0.85, -0.8,
                         -0.75, -0.7, -0.6, -0.3,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0])/(24.0*3600.0)
    # Advective qt forcing  for theta [g/kg/d] converted to [kg/kg/sec]
    r_tend_in = off_arr([ 0.0,   1.2,  2.0,  2.3,   2.2,   2.1,    1.9,   1.7,   1.5,  1.35,  1.22,  1.08,  0.95,
                           0.82,  0.7,  0.6,   0.5,  0.4,   0.3,  0.2,   0.1,  0.05, 0.0025, 0.0012, 0.0006,  0.0,
                           0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0])/(24.0*3600.0)/1000.0
    # Radiative T forcing [K/d] converted to [K/sec]
    Ttend_in = off_arr([ 0.0,  -1.0, -2.2, -3.0,  -3.5,  -3.8,   -4.0,  -4.1,  -4.2,  -4.2,  -4.1,  -4.0, -3.85,
                          -3.7, -3.5, -3.25, -3.0, -2.8,  -2.5, -2.1,  -1.7,  -1.3,   -1.0,   -0.7,   -0.5, -0.4,
                          -0.3, -0.2, -0.1,-0.05,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0])/(24.0*3600.0)
    #! format: on

    Qtend_in = r_tend_in ./ (1 .+ r_tend_in) # convert mixing ratio to specific humidity

    self.Fo.dqtdt = pyinterp(Gr.z_half, z_in, Qtend_in)
    self.Fo.dTdt = pyinterp(Gr.z_half, z_in, Ttend_in) + pyinterp(Gr.z_half, z_in, RAD_in)
end

function initialize_radiation(self::CasesBase{GATE_III}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end

function TurbulenceConvection.initialize_io(self::CasesBase{GATE_III}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{GATE_III}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end

function update_surface(self::CasesBase{GATE_III}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV) # here lhf and shf are needed for calcualtion of bflux in surface and thus u_star
end

function update_forcing(self::CasesBase{GATE_III}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{GATE_III}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

"""
see Stevens et al 2005
Evaluation of Large-Eddy Simulations via Observations of Nocturnal Marine Stratocumulus.
Mon. Wea. Rev., 133, 1443â€“1462.
doi: http://dx.doi.org/10.1175/MWR2930.1
"""
function DYCOMS_RF01(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "DYCOMS_RF01"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceFixedFlux; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingDYCOMS_RF01}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationDYCOMS_RF01}(; Gr, Ref)
    inversion_option = "thetal_maxgrad"
    return TurbulenceConvection.CasesBase{DYCOMS_RF01}(; casename = "DYCOMS_RF01", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{DYCOMS_RF01}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1017.8 * 100.0
    Ref.qtg = 9.0 / 1000.0
    # Use an exner function with values for Rd, and cp given in Stevens 2005 to compute temperature
    Ref.Tg = 289.0 * exner_c(Ref.Pg; kappa = dycoms_Rd / dycoms_cp)
    initialize(Ref, Gr, Stats)
end
# helper function
"""
Compute thetal using constants from Stevens et al 2005 DYCOMS case.
:param p: pressure [Pa]
:param T: temperature [K]
:param ql: liquid water specific humidity
:return: theta l
"""
function dycoms_compute_thetal(self::CasesBase{DYCOMS_RF01}, p_, T_, ql_)
    theta_ = T_ / exner_c(p_, kappa = dycoms_Rd / dycoms_cp)
    return theta_ * exp(-1.0 * dycoms_L * ql_ / (dycoms_cp * T_))
end
# helper function
"""
Use saturation adjustment scheme to compute temperature and ql given thetal and qt.
We can"t use the default TurbulenceConvection function because of different values of cp, Rd and L
:param p: pressure [Pa]
:param thetal: liquid water potential temperature  [K]
:param qt:  total water specific humidity
:return: T, ql
"""
function dycoms_sat_adjst(self::CasesBase{DYCOMS_RF01}, p_, thetal_, qt_)
    #Compute temperature
    t_1 = thetal_ * exner_c(p_, kappa = dycoms_Rd / dycoms_cp)
    #Compute saturation vapor pressure
    pv_star_1 = pv_star(t_1)
    #Compute saturation specific humidity
    qs_1 = qv_star_c(p_, qt_, pv_star_1)

    if qt_ <= qs_1
        #If not saturated return temperature and ql = 0.0
        return t_1, 0.0
    else
        ql_1 = qt_ - qs_1
        f_1 = thetal_ - dycoms_compute_thetal(self, p_, t_1, ql_1)
        t_2 = t_1 + dycoms_L * ql_1 / dycoms_cp
        pv_star_2 = pv_star(t_2)
        qs_2 = qv_star_c(p_, qt_, pv_star_2)
        ql_2 = qt_ - qs_2

        while fabs(t_2 - t_1) >= 1e-9
            pv_star_2 = pv_star(t_2)
            qs_2 = qv_star_c(p_, qt_, pv_star_2)
            ql_2 = qt_ - qs_2
            f_2 = thetal_ - dycoms_compute_thetal(self, p_, t_2, ql_2)
            t_n = t_2 - f_2 * (t_2 - t_1) / (f_2 - f_1)
            t_1 = t_2
            t_2 = t_n
            f_1 = f_2
        end

        return t_2, ql_2
    end
end

function initialize_profiles(self::CasesBase{DYCOMS_RF01}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    thetal = TurbulenceConvection.pyzeros(Gr.nzg) # helper variable to recalculate temperature
    ql = TurbulenceConvection.pyzeros(Gr.nzg) # DYCOMS case is saturated
    qi = 0.0                                             # no ice

    @inbounds for k in real_center_indicies(Gr)
        # thetal profile as defined in DYCOMS
        if Gr.z_half[k] <= 840.0
            thetal[k] = 289.0
        end
        if Gr.z_half[k] > 840.0
            thetal[k] = (297.5 + (Gr.z_half[k] - 840.0)^(1.0 / 3.0))
        end

        # qt profile as defined in DYCOMS
        if Gr.z_half[k] <= 840.0
            GMV.QT.values[k] = 9.0 / 1000.0
        end
        if Gr.z_half[k] > 840.0
            GMV.QT.values[k] = 1.5 / 1000.0
        end

        # ql and T profile
        # (calculated by saturation adjustment using thetal and qt values provided in DYCOMS
        # and using Rd, cp and L constants as defined in DYCOMS)
        GMV.T.values[k], GMV.QL.values[k] = dycoms_sat_adjst(self, Ref.p0_half[k], thetal[k], GMV.QT.values[k])

        # thermodynamic variable profile (either entropy or thetal)
        # (calculated based on T and ql profiles.
        # Here we use Rd, cp and L constants as defined in TurbulenceConvection)
        GMV.THL.values[k] = t_to_thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], GMV.QL.values[k], qi)
        GMV.H.values[k] = t_to_thetali_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], GMV.QL.values[k], qi)


        # buoyancy profile
        qv = GMV.QT.values[k] - qi - GMV.QL.values[k]
        rho = rho_c(Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], qv)
        GMV.B.values[k] = buoyancy_c(Ref.rho0_half[k], rho)

        # velocity profile (geostrophic)
        GMV.U.values[k] = 7.0
        GMV.V.values[k] = -5.5
    end

    # fill out boundary conditions
    set_bcs(GMV.U, Gr)
    set_bcs(GMV.V, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.QL, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.THL, Gr)
    set_bcs(GMV.T, Gr)
    set_bcs(GMV.B, Gr)
    return
end

function initialize_surface(self::CasesBase{DYCOMS_RF01}, Gr::Grid, Ref::ReferenceState)
    self.Sur.zrough = 1.0e-4
    self.Sur.ustar_fixed = false
    self.Sur.cm = 0.0011

    # sensible heat flux
    self.Sur.shf = 15.0
    # latent heat flux
    self.Sur.lhf = 115.0

    self.Sur.Tsurface = 292.5    # K      # i.e. the SST from DYCOMS setup
    self.Sur.qsurface = 13.84e-3 # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?
    #density_surface  = 1.22     # kg/m^3

    # buoyancy flux
    theta_flux = self.Sur.shf / cpm_c(self.Sur.qsurface) / Ref.rho0[Gr.gw - 1]
    qt_flux = self.Sur.lhf / latent_heat(self.Sur.Tsurface) / Ref.rho0[Gr.gw - 1]
    theta_surface = self.Sur.Tsurface / exner_c(Ref.Pg)
    self.Sur.bflux =
        g * (
            (theta_flux + (eps_vi - 1.0) * (theta_surface * qt_flux + self.Sur.qsurface * theta_flux)) /
            (theta_surface * (1.0 + (eps_vi - 1) * self.Sur.qsurface))
        )
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    initialize(self.Sur)
end
function initialize_forcing(self::CasesBase{DYCOMS_RF01}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)

    # geostrophic velocity profiles
    self.Fo.ug .= 7.0
    self.Fo.vg .= -5.5

    # large scale subsidence
    divergence = 3.75e-6    # divergence is defined twice: here and in __init__ of ForcingDYCOMS_RF01 class
    # To be able to have self.Fo.divergence available here,
    # we would have to change the signature of ForcingBase class
    @inbounds for k in real_center_indicies(Gr)
        self.Fo.subsidence[k] = -Gr.z_half[k] * divergence
    end

    # no large-scale drying
    self.Fo.dqtdt .= 0.0 #kg/(kg * s)
end

function initialize_radiation(self::CasesBase{DYCOMS_RF01}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Rad.Gr = Gr
    self.Rad.Ref = Ref
    initialize(self.Rad, GMV)

    # no large-scale drying
    self.Rad.dqtdt .= 0.0 #kg/(kg * s)

    # Radiation based on eq. 3 in Stevens et. al., (2005)
    # cloud-top cooling + cloud-base warming + cooling in free troposphere
    TurbulenceConvection.calculate_radiation(self.Rad, GMV)
end

function TurbulenceConvection.initialize_io(self::CasesBase{DYCOMS_RF01}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
    initialize_io(self.Fo, Stats)
end
function TurbulenceConvection.io(self::CasesBase{DYCOMS_RF01}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
    io(self.Fo, Stats)
end
function update_surface(self::CasesBase{DYCOMS_RF01}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV)
end
function update_forcing(self::CasesBase{DYCOMS_RF01}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end

function update_radiation(self::CasesBase{DYCOMS_RF01}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

function GABLS(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "GABLS"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceMoninObukhovDry; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = true
    latitude = 73.0
    Fo.coriolis_param = 1.39e-4 # s^{-1}
    # Fo.coriolis_param = 2.0 * omega * np.sin(latitude * π / 180.0 ) # s^{-1}
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{GABLS}(; casename = "GABLS", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{GABLS}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1.0e5  #Pressure at ground
    Ref.Tg = 265.0  #Temperature at ground
    Ref.qtg = 1.0e-12 #Total water mixing ratio at TurbulenceConvection. if set to 0, alpha0, rho0, p0 are NaN (TBD)
    initialize(Ref, Gr, Stats)
end
function initialize_profiles(self::CasesBase{GABLS}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    thetal = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0 # IC of GABLS cloud-free

    @inbounds for k in real_center_indicies(Gr)
        #Set wind velocity profile
        GMV.U.values[k] = 8.0
        GMV.V.values[k] = 0.0

        #Set Thetal profile
        if Gr.z_half[k] <= 100.0
            thetal[k] = 265.0
        else
            thetal[k] = 265.0 + (Gr.z_half[k] - 100.0) * 0.01
        end

        #Set qt profile
        GMV.QT.values[k] = 0.0
    end

    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = thetal[k]
        GMV.T.values[k] = thetal[k] * exner_c(Ref.p0_half[k]) # No water content
        GMV.THL.values[k] = thetal[k]
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.V, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{GABLS}, Gr::Grid, Ref::ReferenceState)
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.zrough = 0.1
    self.Sur.Tsurface = 265.0
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{GABLS}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    @inbounds for k in real_center_indicies(Gr)
        # Geostrophic velocity profiles.
        self.Fo.ug[k] = 8.0
        self.Fo.vg[k] = 0.0
    end
    return nothing
end

function initialize_radiation(self::CasesBase{GABLS}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end

function TurbulenceConvection.initialize_io(self::CasesBase{GABLS}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end

function TurbulenceConvection.io(self::CasesBase{GABLS}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end

function update_surface(self::CasesBase{GABLS}, GMV::GridMeanVariables, TS::TimeStepping)
    self.Sur.Tsurface = 265.0 - (0.25 / 3600.0) * TS.t
    update(self.Sur, GMV)
end

function update_forcing(self::CasesBase{GABLS}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{GABLS}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

# Not fully implemented yet - Ignacio
function SP(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "SP"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceSullivanPatton; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingStandard}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "critical_Ri"
    Fo.apply_coriolis = true
    Fo.coriolis_param = 1.0e-4 # s^{-1}
    # Fo.coriolis_param = 2.0 * omega * np.sin(latitude * π / 180.0 ) # s^{-1}
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{SP}(; casename = "SP", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{SP}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1.0e5  #Pressure at ground
    Ref.Tg = 300.0  #Temperature at ground
    Ref.qtg = 1.0e-4   #Total water mixing ratio at TurbulenceConvection. if set to 0, alpha0, rho0, p0 are NaN.
    initialize(Ref, Gr, Stats)
end

function initialize_profiles(self::CasesBase{SP}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)
    thetal = TurbulenceConvection.pyzeros(Gr.nzg)
    ql = 0.0
    qi = 0.0 # IC of SP cloud-free

    @inbounds for k in real_center_indicies(Gr)
        GMV.U.values[k] = 1.0
        GMV.V.values[k] = 0.0
        #Set Thetal profile
        if Gr.z_half[k] <= 974.0
            thetal[k] = 300.0
        elseif Gr.z_half[k] < 1074.0
            thetal[k] = 300.0 + (Gr.z_half[k] - 974.0) * 0.08
        else
            thetal[k] = 308.0 + (Gr.z_half[k] - 1074.0) * 0.003
        end

        #Set qt profile
        GMV.QT.values[k] = 0.0
    end

    @inbounds for k in real_center_indicies(Gr)
        GMV.H.values[k] = thetal[k]
        GMV.T.values[k] = thetal[k] * exner_c(Ref.p0_half[k])
        GMV.THL.values[k] = thetal[k]
    end

    set_bcs(GMV.U, Gr)
    set_bcs(GMV.V, Gr)
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    set_bcs(GMV.T, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{SP}, Gr::Grid, Ref::ReferenceState)
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.zrough = 0.1
    self.Sur.Tsurface = 300.0
    theta_surface = self.Sur.Tsurface / exner_c(Ref.Pg)
    theta_flux = 0.24
    self.Sur.bflux = g * theta_flux / theta_surface
    # self.Sur.bflux = 0.24 * exner_c(Ref.p0_half[Gr.gw]) * g / (Ref.p0_half[Gr.gw]*Ref.alpha0_half[Gr.gw]/Rd)
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{SP}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
    @inbounds for k in real_center_indicies(Gr)
        # Geostrophic velocity profiles. vg = 0
        self.Fo.ug[k] = 1.0
        self.Fo.vg[k] = 0.0
    end
end
function initialize_radiation(self::CasesBase{SP}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end
function TurbulenceConvection.initialize_io(self::CasesBase{SP}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end
function TurbulenceConvection.io(self::CasesBase{SP}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end
function update_surface(self::CasesBase{SP}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV)
end
function update_forcing(self::CasesBase{SP}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{SP}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

function DryBubble(namelist, Gr::Grid, Ref::ReferenceState)
    casename = "DryBubble"
    Sur = TurbulenceConvection.SurfaceBase(TurbulenceConvection.SurfaceNone; Gr, Ref, namelist)
    Fo = TurbulenceConvection.ForcingBase{TurbulenceConvection.ForcingNone}(; Gr, Ref)
    Rad = TurbulenceConvection.RadiationBase{TurbulenceConvection.RadiationNone}(; Gr, Ref)
    inversion_option = "theta_rho"
    Fo.apply_coriolis = false
    Fo.apply_subsidence = false
    return TurbulenceConvection.CasesBase{DryBubble}(; casename = "DryBubble", inversion_option, Sur, Fo, Rad)
end

function initialize_reference(self::CasesBase{DryBubble}, Gr::Grid, Ref::ReferenceState, Stats::NetCDFIO_Stats)
    Ref.Pg = 1.0e5  #Pressure at ground
    Ref.qtg = 1.0e-5
    Ref.Tg = 296
    initialize(Ref, Gr, Stats)
end

function initialize_profiles(self::CasesBase{DryBubble}, Gr::Grid, GMV::GridMeanVariables, Ref::ReferenceState)

    @inbounds for k in center_indicies(Gr)
        GMV.U.values[k] = 0.01
    end

    n_updrafts = 1
    # initialize Grid Mean Profiles of thetali and qt
    #! format: off
    z_in = off_arr([
                      25.,   75.,  125.,  175.,  225.,  275.,  325.,  375.,  425.,
                     475.,  525.,  575.,  625.,  675.,  725.,  775.,  825.,  875.,
                     925.,  975., 1025., 1075., 1125., 1175., 1225., 1275., 1325.,
                    1375., 1425., 1475., 1525., 1575., 1625., 1675., 1725., 1775.,
                    1825., 1875., 1925., 1975., 2025., 2075., 2125., 2175., 2225.,
                    2275., 2325., 2375., 2425., 2475., 2525., 2575., 2625., 2675.,
                    2725., 2775., 2825., 2875., 2925., 2975., 3025., 3075., 3125.,
                    3175., 3225., 3275., 3325., 3375., 3425., 3475., 3525., 3575.,
                    3625., 3675., 3725., 3775., 3825., 3875., 3925., 3975., 4025.,
                    4075., 4125., 4175., 4225., 4275., 4325., 4375., 4425., 4475.,
                    4525., 4575., 4625., 4675., 4725., 4775., 4825., 4875., 4925.,
                    4975., 5025., 5075., 5125., 5175., 5225., 5275., 5325., 5375.,
                    5425., 5475., 5525., 5575., 5625., 5675., 5725., 5775., 5825.,
                    5875., 5925., 5975., 6025., 6075., 6125., 6175., 6225., 6275.,
                    6325., 6375., 6425., 6475., 6525., 6575., 6625., 6675., 6725.,
                    6775., 6825., 6875., 6925., 6975., 7025., 7075., 7125., 7175.,
                    7225., 7275., 7325., 7375., 7425., 7475., 7525., 7575., 7625.,
                    7675., 7725., 7775., 7825., 7875., 7925., 7975., 8025., 8075.,
                    8125., 8175., 8225., 8275., 8325., 8375., 8425., 8475., 8525.,
                    8575., 8625., 8675., 8725., 8775., 8825., 8875., 8925., 8975.,
                    9025., 9075., 9125., 9175., 9225., 9275., 9325., 9375., 9425.,
                    9475., 9525., 9575., 9625., 9675., 9725., 9775., 9825., 9875.,
                    9925., 9975.
    ])
    thetali_in = off_arr([
                    299.9834, 299.9836, 299.9841, 299.985 , 299.9864, 299.9883,
                    299.9907, 299.9936, 299.9972, 300.0012, 300.0058, 300.011 ,
                    300.0166, 300.0228, 300.0293, 300.0363, 300.0436, 300.0512,
                    300.0591, 300.0672, 300.0755, 300.0838, 300.0921, 300.1004,
                    300.1086, 300.1167, 300.1245, 300.132 , 300.1393, 300.1461,
                    300.1525, 300.1583, 300.1637, 300.1685, 300.1726, 300.1762,
                    300.179 , 300.1812, 300.1826, 300.1833, 300.1833, 300.1826,
                    300.1812, 300.179 , 300.1762, 300.1727, 300.1685, 300.1637,
                    300.1584, 300.1525, 300.1461, 300.1393, 300.1321, 300.1245,
                    300.1167, 300.1087, 300.1005, 300.0922, 300.0838, 300.0755,
                    300.0673, 300.0592, 300.0513, 300.0437, 300.0364, 300.0294,
                    300.0228, 300.0167, 300.0111, 300.0059, 300.0013, 299.9972,
                    299.9937, 299.9908, 299.9884, 299.9865, 299.9851, 299.9842,
                    299.9837, 299.9835, 299.9835, 299.9835, 299.9835, 299.9835,
                    299.9835, 299.9835, 299.9835, 299.9835, 299.9835, 299.9835,
                    299.9835, 299.9835, 299.9835, 299.9835, 299.9835, 299.9835,
                    299.9835, 299.9835, 299.9835, 299.9835, 299.9835, 299.9835,
                    299.9835, 299.9835, 299.9835, 299.9835, 299.9835, 299.9835,
                    299.9835, 299.9835, 299.9835, 299.9835, 299.9835, 299.9835,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9836, 299.9836,
                    299.9836, 299.9836, 299.9836, 299.9836, 299.9837, 299.9837,
                    299.9837, 299.9837, 299.9837, 299.9837, 299.9837, 299.9837,
                    299.9837, 299.9837, 299.9837, 299.9837, 299.9837, 299.9837,
                    299.9837, 299.9837, 299.9837, 299.9837, 299.9837, 299.9837,
                    299.9837, 299.9837
    ])
    #! format: on
                       #LES temperature_mean in K
    thetali = TurbulenceConvection.pyzeros(Gr.nzg)
    z_half_in = off_arr(Gr.z_half[(Gr.gw):(Gr.nzg - Gr.gw)])
    thetali[(Gr.gw):(Gr.nzg - Gr.gw)] = pyinterp(z_half_in, z_in, thetali_in)
    GMV.THL.values .= thetali
    GMV.H.values .= thetali
    @inbounds for k in real_center_indicies(Gr)
        GMV.QT.values[k] = 0.0
    end
    set_bcs(GMV.QT, Gr)
    set_bcs(GMV.H, Gr)
    satadjust(GMV)
end

function initialize_surface(self::CasesBase{DryBubble}, Gr::Grid, Ref::ReferenceState)
    self.Sur.Gr = Gr
    self.Sur.Ref = Ref
    self.Sur.qsurface = 1.0e-5
    self.Sur.shf = 8.0e-3 * cpm_c(self.Sur.qsurface) * Ref.rho0[Gr.gw - 1]
    initialize(self.Sur)
end

function initialize_forcing(self::CasesBase{DryBubble}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    self.Fo.Gr = Gr
    self.Fo.Ref = Ref
    initialize(self.Fo, GMV)
end

function initialize_radiation(self::CasesBase{DryBubble}, Gr::Grid, Ref::ReferenceState, GMV::GridMeanVariables)
    initialize(self.Rad, GMV)
end

function TurbulenceConvection.initialize_io(self::CasesBase{DryBubble}, Stats::NetCDFIO_Stats)
    initialize_io(self, Stats, BaseCase())
end

function TurbulenceConvection.io(self::CasesBase{DryBubble}, Stats::NetCDFIO_Stats)
    io(self, Stats, BaseCase())
end

function update_surface(self::CasesBase{DryBubble}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Sur, GMV)
end

function update_forcing(self::CasesBase{DryBubble}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Fo, GMV)
end
function update_radiation(self::CasesBase{DryBubble}, GMV::GridMeanVariables, TS::TimeStepping)
    update(self.Rad, GMV)
end

end
