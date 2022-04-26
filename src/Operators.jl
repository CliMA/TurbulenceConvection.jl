#####
##### Masked operators
#####

"""
    shrink_mask!(f, mask)

Shrinks the subdomains of `f` where `mask == 1`.

Example:

```julia
using Test
f = Int[0, 0, 0, 1, 1, 1, 0, 0, 1, 1]
mask = Bool[0, 0, 0, 1, 1, 1, 0, 0, 1, 1]
@test shrink_mask!(f, mask) == Bool[0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
```
"""
function shrink_mask!(f::AbstractArray, mask::AbstractArray)
    @inbounds for (i, m) in enumerate(mask)
        if i == 1 || i == length(mask)
            f[i] = m
        elseif m == 1 && mask[i - 1] == 0
            f[i] = 0
        elseif m == 1 && mask[i + 1] == 0
            f[i] = 0
        else
            f[i] = m
        end
    end
end

function special_interp(var, zc, z_star::FT) where {FT}
    LBF = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))
    zero_bcs = (; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    I0f = CCO.InterpolateC2F(; zero_bcs...)
    is_updraft_top = @. Int(zc.z ≈ z_star) * FT(1)
    return @. LBF(var*is_updraft_top) + I0f(var*(1-is_updraft_top))
end
