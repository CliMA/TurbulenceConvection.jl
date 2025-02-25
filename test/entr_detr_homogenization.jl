
FT = Float64

N_iter = 100

x = FT[1, 2]
x = FT[2, 1]

x0 = copy(x)

c_entr = 0.1
c_detr = 0.5

x_out = zeros(FT, (2, N_iter+1))
x_out[:, 1] .= x


@info "Iter: 0, x: $x"
Δt = 0.5
for i in 1:N_iter
    x .+= [ -c_detr * x[1] * Δt + c_entr * x[2] * Δt, c_detr * x[1] * Δt - c_entr * x[2] * Δt]
    @info "Iter: $i, x: $x"
    
    x_out[:, i+1] .= x
end

k_1 = (x0[1] + x0[2]) / (c_entr + c_detr)
k_2 = (x0[1]*c_detr - x0[2]*c_entr) / (c_entr + c_detr)
x_1_analytical(t) = k_1 * c_entr + k_2 * exp.(-(c_entr + c_detr)*t)
x_2_analytical(t) = k_1 * c_detr - k_2 * exp.(-(c_entr + c_detr)*t)
t = 0:Δt:N_iter*Δt

x_0_steady = [k_1 * c_entr, k_1 * c_detr]

using Plots
thisdir = dirname(@__FILE__)
plot(; size=(800, 400), dpi=600, legend=:left)
plot!(t, x_out[1, :], color=:red, label="x1")
plot!(t, x_out[2, :], color=:blue, label="x2")
plot!(t, x_1_analytical.(t), color=:red, linestyle=:dash, label="x1 analytical")
plot!(t, x_2_analytical.(t), color=:blue, linestyle=:dash, label="x2 analytical")
hline!([x_0_steady[1]], color=:red, linestyle=:dot, label="x1 steady state")
hline!([x_0_steady[2]], color=:blue, linestyle=:dot, label="x2 steady state")
plot!(t, x_out[1,:] + x_out[2,:], color=:black, linestyle=:dash, label="x1 + x2")
savefig(joinpath(thisdir, "entr_detr_homogenization.png"))




plot(; size=(800, 400), dpi=600, legend=:left)
plot!(1:N_iter+1, x_out[1, :] ./ x_out[2, :], label="ratio")
plot!(1:N_iter+1, x_out[1,1]/x_out[2,1] * exp.(-(0:N_iter)) * (1 - x_out[1,1]/x_out[2,1]), color=:red, linestyle=:dot, label="x2/x1 * exp(-N_iter)")
hline!([c_entr / c_detr], label="c_entr / c_detr", color=:black, linestyle=:dash)
savefig(joinpath(thisdir, "entr_detr_homogenization_ratio.png"))