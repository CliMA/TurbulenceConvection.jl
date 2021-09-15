# TODO: check that these primitives are correct

# Use stop-1 to emulate python ranges
xrange(start, stop, step = 1) = range(start + 1, stop; step = step)
xrange(stop) = xrange(0, stop)

off_arr(a::AbstractArray) = a

function pyinterp(x, xp, fp)
    spl = Dierckx.Spline1D(xp, fp; k = 1)
    return spl(vec(x))
end
