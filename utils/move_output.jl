if haskey(ENV, "BUILDKITE_COMMIT") && haskey(ENV, "BUILDKITE_BRANCH")
    commit = ENV["BUILDKITE_COMMIT"]
    branch = ENV["BUILDKITE_BRANCH"]
    # Note: cluster_data_prefix is also defined in integration_tests/utils/compute_mse.jl
    cluster_data_prefix = "/central/scratch/climaci/turbulenceconvection-main"

    @info "pwd() = $(pwd())"
    @info "branch = $(branch)"
    @info "commit = $(commit)"

    using Glob
    if branch == "staging"
        commit_sha = commit[1:7]
        path = joinpath(cluster_data_prefix, commit_sha)
        mkpath(path)
        for folder_name in glob("Output.*")
            src = folder_name
            dst = joinpath(path, folder_name)
            @info "Moving $src to $dst"
            mv(src, dst; force = true)
        end
        @info "readdir(): $(readdir(path))"
    end
else
    @info "ENV keys: $(keys(ENV))"
end
