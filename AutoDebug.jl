try
    include("RunDycoms.jl")
catch
    include("DebugDycoms.jl")
end
