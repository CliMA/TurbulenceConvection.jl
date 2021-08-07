using OffsetArrays

# TODO: check that these primitives are correct
fmax(a, b) = max(a, b)
fmin(a, b) = min(a, b)
fabs(a) = abs(a)

argwhere(a) = findall(a .> 0)
# Use stop-1 to emulate python ranges
xrange(start, stop, step = 1) = range(start + 1, stop; step = step)
xrange(stop) = xrange(0, stop)

off_arr(a::AbstractArray) = a

pow(a::Real, b::Real) = a^b
pow(a::AbstractVector, b::AbstractVector) = a .^ b
# TODO: double check this translation (from numpy)
power(a::Real, b::Real) = a^b
power(a::AbstractVector, b::AbstractVector) = a .^ b

function pyinterp(x::T, xp::T, fp::T) where {T}
    spl = Dierckx.Spline1D([xp...], [fp...]; k = 1)
    return spl([x...])
end
