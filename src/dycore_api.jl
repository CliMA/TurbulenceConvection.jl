#####
##### Dycore API
#####

abstract type FieldLocation end
struct CentField <: FieldLocation end
struct FaceField <: FieldLocation end

field_loc(::CentField) = :cent
field_loc(::FaceField) = :face

#=
This file provides a list of methods that TurbulenceConvection.jl
expects that the host dycore supports. This is experimental, as
we're not sure how the data structures / flow control will shake out.
=#

#####
##### Grid mean fields
#####

""" Prognostic fields for the host model """
prognostic(state, fl) = getproperty(state.prog, field_loc(fl))
center_prog_grid_mean(state) = prognostic(state, CentField())
face_prog_grid_mean(state) = prognostic(state, FaceField())

""" Auxiliary fields for the host model """
aux(state, fl) = getproperty(state.aux, field_loc(fl))
center_aux_grid_mean(state) = aux(state, CentField())

""" Tendency fields for the host model """
tendencies(state, fl) = getproperty(state.tendencies, field_loc(fl))
center_tendencies_grid_mean(state) = tendencies(state, CentField())

""" Reference state fields for the host model """
ref_state(state, fl) = aux(state, fl).ref_state
face_ref_state(state) = ref_state(state, FaceField())
center_ref_state(state) = ref_state(state, CentField())

#####
##### TurbulenceConvection fields
#####

#= Prognostic fields for TurbulenceConvection =#
prognostic_tc(state, fl) = prognostic(state, fl).turbconv
center_prog_updrafts(state) = prognostic_tc(state, CentField()).up
face_prog_updrafts(state) = prognostic_tc(state, FaceField()).up
center_prog_environment(state) = prognostic_tc(state, CentField()).en
face_prog_environment(state) = prognostic_tc(state, FaceField()).en

#= Auxiliary fields for TurbulenceConvection =#
aux_turbconv(state, fl) = aux(state, fl).turbconv
center_aux_tc(state) = aux_turbconv(state, CentField())
face_aux_tc(state) = aux_turbconv(state, FaceField())
center_aux_updrafts(state) = aux_turbconv(state, CentField()).up
face_aux_updrafts(state) = aux_turbconv(state, FaceField()).up
center_aux_environment(state) = aux_turbconv(state, CentField()).en
face_aux_environment(state) = aux_turbconv(state, FaceField()).en

#= Tendency fields for TurbulenceConvection =#
tendencies_tc(state, fl) = tendencies(state, fl).turbconv
center_tendencies_tc(state) = tendencies_tc(state, CentField())
face_tendencies_tc(state) = tendencies_tc(state, FaceField())
center_tendencies_updrafts(state) = tendencies_tc(state, CentField()).up
face_tendencies_updrafts(state) = tendencies_tc(state, FaceField()).up
