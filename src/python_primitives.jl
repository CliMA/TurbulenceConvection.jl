using OffsetArrays

# TODO: check that these primitives are correct
fmax(a, b) = max(a, b)
fmin(a, b) = min(a, b)
fabs(a) = abs(a)

argwhere(a) = findall(a .> 0)
# Use stop-1 to emulate python ranges
xrange(start, stop, step = 1) = range(start, stop - 1; step = step)
xrange(stop) = xrange(0, stop)

revxrange(start, stop, step = 1) = range(start, stop; step = step)
revxrange(stop) = revxrange(0, stop)

pyzeros(n::Int) = OffsetArray(zeros(n), 0:(n - 1))
pyzeros(m::Int, n::Int) = OffsetArray(zeros(m, n), 0:(m - 1), 0:(n - 1))

pyones(n::Int) = OffsetArray(ones(n), 0:(n - 1))
pyones(m::Int, n::Int) = OffsetArray(ones(m, n), 0:(m - 1), 0:(n - 1))

function off_arr(a::AbstractArray)
    dims = ntuple(ndims(a)) do i
        0:(size(a, i) - 1)
    end
    return OffsetArray(a, dims)
end

pow(a::Real, b::Real) = a^b
pow(a::AbstractVector, b::AbstractVector) = a .^ b
# TODO: double check this translation (from numpy)
power(a::Real, b::Real) = a^b
power(a::AbstractVector, b::AbstractVector) = a .^ b
linspace(a, b; num = 50) = range(a, b; length = num)

function pyinterp(x::T, xp::T, fp::T) where {T}
    spl = Dierckx.Spline1D([xp...], [fp...]; k = 1)
    return off_arr(spl([x...]))
end
