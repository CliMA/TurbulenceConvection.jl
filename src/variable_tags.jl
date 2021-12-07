abstract type AbstractVariableTag end
abstract type AbstractBulkVariableTag <: AbstractVariableTag end

struct BulkEntr <: AbstractBulkVariableTag end
struct BulkDetr <: AbstractBulkVariableTag end
struct BulkAspRatio <: AbstractBulkVariableTag end
struct BulkFracTurbEntr <: AbstractBulkVariableTag end
struct BulkHorizKEddy <: AbstractBulkVariableTag end

tagname(::BulkEntr) = "entrainment_sc"
tagname(::BulkDetr) = "detrainment_sc"
tagname(::BulkAspRatio) = "asp_ratio"
tagname(::BulkFracTurbEntr) = "turbulent_entrainment"
tagname(::BulkHorizKEddy) = "horiz_K_eddy"

tagdims(::AbstractVariableTag) = ("zc", "t")
taggroup(::AbstractVariableTag) = "profiles"

function initialize_io_new(Stats::NetCDFIO_Stats)
    for tag in (
            BulkEntr(),
            BulkDetr(),
            BulkAspRatio(),
            BulkFracTurbEntr(),
            BulkHorizKEddy(),
        )
        add_field(Stats, tagname(tag); dims = tagdims(tag), group = taggroup(tag))
    end
end

#####
##### Generic methods
#####

function compute_export!(tag::AbstractVariableTag, state::State, grid, edmf, Stats)
    term = center_aux_turbconv(state).ϕ_temporary
    compute!(tag, term, state, grid, edmf)
    write_field(Stats, tagname(tag), term; group = taggroup(tag))
end

function compute!(tag::AbstractVariableTag, term::CC.Fields.FiniteDifferenceField, state::State, grid, edmf)
    parent(term) .= 0
    🌐 = (;
        aux_up = center_aux_updrafts(state),
        aux_en = center_aux_environment(state),
        aux_bulk = center_aux_bulk(state),
    )
    compute!(tag, term, 🌐, edmf)
end

function compute!(tag::AbstractBulkVariableTag, term::CC.Fields.FiniteDifferenceField, 🌐, edmf)
    @inbounds for i in 1:n_updrafts(edmf)
        compute∑!(tag, term, 🌐, i)
    end
end

is_updraft(a_bulk) = a_bulk > 0
ife0(B, x) = ifelse(B, x, 0)
compute∑!(tag::AbstractBulkVariableTag, term::CC.Fields.FiniteDifferenceField, 🌐, i::Int) =
    (term .+= ife0.(is_updraft.(🌐.aux_bulk.area), field(tag, 🌐, i)))

#####
##### Individual expressions
#####

field(::BulkEntr, 🌐, i::Int) = @. 🌐.aux_up[i].area * 🌐.aux_up[i].entr_sc / 🌐.aux_bulk.area
field(::BulkDetr, 🌐, i::Int) = @. 🌐.aux_up[i].area * 🌐.aux_up[i].detr_sc / 🌐.aux_bulk.area
field(::BulkAspRatio, 🌐, i::Int) = @. 🌐.aux_up[i].area * 🌐.aux_up[i].asp_ratio / 🌐.aux_bulk.area
field(::BulkFracTurbEntr, 🌐, i::Int) = @. 🌐.aux_up[i].area * 🌐.aux_up[i].frac_turb_entr / 🌐.aux_bulk.area
field(::BulkHorizKEddy, 🌐, i::Int) = @. 🌐.aux_up[i].area * 🌐.aux_up[i].horiz_K_eddy / 🌐.aux_bulk.area


