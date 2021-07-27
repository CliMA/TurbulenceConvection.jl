try
    include("RunCase.jl")
catch
    include("DebugCase.jl")
end
