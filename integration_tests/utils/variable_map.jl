#! format: off
var_map(s::String) = var_map(Val(Symbol(s)))
var_map(s::Symbol) = var_map(Val(s))
var_map(::Val{T}) where {T} = nothing

var_map(::Val{:rho}) = ("rho", ())
var_map(::Val{:u_mean}) = ("u_mean", (:ρ,))
var_map(::Val{:v_mean}) = ("v_mean", (:ρ,))
var_map(::Val{:qt_mean}) = ("qt_mean", (:ρ,))
var_map(::Val{:updraft_fraction}) = ("updraft_fraction", (:ρ,))
var_map(::Val{:updraft_w}) = ("updraft_w", (:ρ, :a))
var_map(::Val{:updraft_qt}) = ("updraft_qt", (:ρ, :a))
var_map(::Val{:updraft_thetali}) = ("updraft_thetali", (:ρ, :a))
var_map(::Val{:tke_mean}) = ("tke_mean", (:ρ, :a))
var_map(::Val{:env_thetali2}) = ("env_thetali2", (:ρ, :a))
var_map(::Val{:env_qt2}) = ("env_qt2", (:ρ, :a))
#! format: on

var_map(::Val{:q_tot_gm}) = ("qt_mean", ())
var_map(::Val{:a_up}) = ("updraft_fraction", ())
var_map(::Val{:w_up}) = ("updraft_w", ())
var_map(::Val{:q_tot_up}) = ("updraft_qt", ())
var_map(::Val{:θ_liq_up}) = ("updraft_thetali", ())
var_map(::Val{:tke_en}) = ("tke_mean", ())
#! format: on
