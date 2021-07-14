#! format: off
var_map(s::String) = var_map(Val(Symbol(s)))
var_map(s::Symbol) = var_map(Val(s))
var_map(::Val{T}) where {T} = nothing

# var_map(::Val{:tc_var}) = "les_var"
var_map(::Val{:rho}) = "rho"
var_map(::Val{:u_mean}) = "u_mean"
var_map(::Val{:v_mean}) = "v_mean"
var_map(::Val{:qt_mean}) = "qt_mean"
# var_map(::Val{:updraft_fraction}) = "updraft_area"
var_map(::Val{:updraft_area}) = "updraft_fraction"
var_map(::Val{:updraft_w}) = "updraft_w"
var_map(::Val{:updraft_qt}) = "updraft_qt"
var_map(::Val{:updraft_ql}) = "updraft_ql"
var_map(::Val{:updraft_qr}) = "updraft_ql"
var_map(::Val{:updraft_thetal}) = "updraft_thetali"
var_map(::Val{:tke_mean}) = "tke_mean"
#! format: on
