using DataFrames

function compute_chimeric(df::DataFrame, chimeric_col::Symbol, sb_col::Symbol, umi_col::Symbol)::Nothing
    if nrow(df) == 0
        df[!, chimeric_col] = Bool[]
        return nothing
    end

    print_start("Computing $(sb_col) chimeras... ")

    sort!(df, [sb_col, umi_col, :reads], rev = [false, false, true])
    before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [sb_col, umi_col]]))
    after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [sb_col, umi_col, :reads]]), false)

    df[!, chimeric_col] = before_same .| after_same

    println_done()
    return nothing
end

function compute_chimeric1(df::DataFrame)::Nothing
    compute_chimeric(df, :chimeric1, :sb1_i, :umi1_i)
    return nothing
end

function compute_chimeric2(df::DataFrame)::Nothing
    compute_chimeric(df, :chimeric2, :sb2_i, :umi2_i)
    return nothing
end

function remove_chimeras(df::DataFrame, metadata::Dict{String,Int64})::Nothing
    print_start("Removing chimeras... ")
    metadata["R1_chimeric"] = sum(df.chimeric1)
    metadata["R2_chimeric"] = sum(df.chimeric2)

    subset!(df, :chimeric1 => x -> .!x, :chimeric2 => x -> .!x)
    select!(df, Not([:chimeric1, :chimeric2]))
    metadata["umis_chimeric"] = metadata["umis_exact"] - nrow(df)

    println_done()
    return nothing
end

function write_chimeric(df::DataFrame, out_path::String, col::Symbol, file_suffix::String = "")::Nothing
    df_chimeric = df[df[!, col], [:sb1_i, :umi1_i, :sb2_i, :umi2_i]]
    write_df(df_chimeric, joinpath(out_path, "$(col)$(file_suffix).csv.gz"))
    return nothing
end

function write_chimeric1(df::DataFrame, out_path::String, file_suffix::String = "")::Nothing
    write_chimeric(df, out_path, :chimeric1, file_suffix)
    return nothing
end

function write_chimeric2(df::DataFrame, out_path::String, file_suffix::String = "")::Nothing
    write_chimeric(df, out_path, :chimeric2, file_suffix)
    return nothing
end

function find_chimerics(in_path::String, col::Symbol)::Vector{String}
    # Load the chimeric paths
    chimerics = readdir(in_path, join=true)
    chimerics = filter(chimeric -> startswith(basename(chimeric), "$col."), chimerics)
    chimerics = filter(chimeric -> endswith(chimeric, ".csv.gz"), chimerics)
    println("$(col)s: ", basename.(chimerics))
    @assert length(chimerics) > 0 "ERROR: No $col files found"
    println("$(length(chimerics)) $col file(s) found\n")
    return chimerics
end

function read_chimeric(df::DataFrame, in_path::String, chimeric_col::Symbol)::Nothing
    chimerics = find_chimerics(in_path, chimeric_col)
    df_chimeric = read_dfs(chimerics, Dict(:sb1_i => UInt64, :umi1_i => UInt64, :sb2_i => UInt64, :umi2_i => UInt64))
    unique!(df_chimeric)
    df_chimeric[!, chimeric_col] = trues(nrow(df_chimeric))
    leftjoin!(df, df_chimeric, on=[:sb1_i, :umi1_i, :sb2_i, :umi2_i])
    df[!, chimeric_col] = coalesce.(df[!, chimeric_col], false)
    return nothing
end

function read_chimerics(df::DataFrame, in_path::String)::Nothing
    read_chimeric(df, in_path, :chimeric1)
    read_chimeric(df, in_path, :chimeric2)
    return nothing
end
