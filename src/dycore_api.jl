#####
##### Dycore API
#####

#=
This file provides a list of methods that TurbulenceConvection.jl
expects that the host dycore supports. This is experimental, as
we're not sure how the data structures / flow control will shake out.
=#

""" The cell center fields of any state vector """
center_prog(state) = state.prog.cent

""" The cell face fields of any state vector """
face_prog(state) = state.prog.face

""" The cell center fields of the edmf updrafts """
center_prog_updrafts(state) = center_prog(state).turbconv.up

""" The cell face fields of the edmf updrafts """
face_prog_updrafts(state) = face_prog(state).turbconv.up

""" The cell center fields of the edmf environment """
center_prog_environment(state) = center_prog(state).turbconv.up

""" The cell face fields of the edmf environment """
face_prog_environment(state) = face_prog(state).turbconv.up

""" The cell center reference state fields """
center_ref_state(state) = state.aux.cent.ref_state

""" The cell face reference state fields """
face_ref_state(state) = state.aux.face.ref_state
