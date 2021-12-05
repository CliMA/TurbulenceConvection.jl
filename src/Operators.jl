#####
##### Masked operators
#####

"""
    shrink_mask(mask)

Shrinks the subdomains where `mask == 1`.

Example:

```julia
using Test
mask = Bool[0, 0, 0, 1, 1, 1, 0, 0, 1, 1]
shrunken_mask = Bool[0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
@test shrink_mask(mask) == shrunken_mask
```
"""
function shrink_mask(mask)
    return map(enumerate(mask)) do (i, m)
        if i == 1 || i == length(mask)
            m
        elseif m == 1 && mask[i - 1] == 0
            m = 0
        elseif m == 1 && mask[i + 1] == 0
            m = 0
        else
            m
        end
    end
end
