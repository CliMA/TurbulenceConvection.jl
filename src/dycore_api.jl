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

""" The cell center fields of any state vector """
prognostic(state, fl) = getproperty(state.prog, field_loc(fl))

""" The cell center reference state fields """
aux(state, fl) = getproperty(state.aux, field_loc(fl))

""" The cell center reference state fields """
ref_state(state, fl) = aux(state, fl).ref_state
face_ref_state(state) = ref_state(state, FaceField())
center_ref_state(state) = ref_state(state, CentField())

#####
##### TC fields
#####

""" The cell center fields of the edmf updrafts """
prognostic_tc(state, fl) = prognostic(state, fl).turbconv
center_prog_updrafts(state) = prognostic_tc(state, CentField()).up
face_prog_updrafts(state) = prognostic_tc(state, FaceField()).up
center_prog_environment(state) = prognostic_tc(state, CentField()).en
face_prog_environment(state) = prognostic_tc(state, FaceField()).en
