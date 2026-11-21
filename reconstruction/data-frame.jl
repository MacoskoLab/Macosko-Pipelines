using CSV
using CodecZlib
using DataFrames

# A df contains full information for pairs of barcodes.
# Each file possibly includes the number of reads or umi.
# Much larger than a data table.

# CSV.jl's compress=true gzips at the default level, which made writing the big frames take
# longer than computing them; level 1 is several times faster for a slightly larger file
function write_csv(path::String, df::DataFrame; compress::Bool = endswith(path, ".gz"))::Nothing
    if compress
        open(path, "w") do file
            io = GzipCompressorStream(file; level=1)
            CSV.write(io, df, writeheader=true)
            close(io)
        end
    else
        CSV.write(path, df, writeheader=true)
    end
    return nothing
end

function write_df(df::DataFrame, path::String; compress::Union{Bool,Nothing} = nothing)::Nothing
    print_start("Writing $(basename(path))... ")
    if compress === nothing
        compress = endswith(path, ".gz")
    end
    write_csv(path, df; compress=compress)
    println_done()
    return nothing
end

function read_df(path::String, types::Dict{Symbol,DataType})::DataFrame
    print_start("Reading $(basename(path))... ")
    df = CSV.read(path, DataFrame; header=true, types=types)
    println_done()
    return df
end

function sort_df(df::DataFrame, cols)::Nothing
    print_start("Sorting by $cols... ")
    sort!(df, cols)
    println_done()
    return nothing
end

function read_dfs(paths::Vector{String}, types::Dict{Symbol,DataType})::DataFrame
    df = reduce(
        (df1, path) -> vcat(df1, read_df(path, types)),
        paths[2:end];
        init = read_df(paths[1], types),
    )
    return df
end

# A tab is only a summary data table with counts of rows from a data frame.
# Much smaller than a data frame.

function write_tab(
    df::DataFrame,
    path::String,
    input_col::Symbol,    # Name of the column in df to count (e.g., :reads)
    key_col::Symbol,      # Name for the key column in output (e.g., :reads_per_umi)
    value_col::Symbol;    # Name for the value column in output (e.g., :umis)
    sort_by::Union{Symbol,Nothing} = nothing,
    compress::Union{Bool,Nothing} = nothing,
)::Nothing
    tab_dict = countmap(df[!, input_col])
    tab_df = DataFrame(
        key_col => collect(keys(tab_dict)),
        value_col => collect(values(tab_dict)),
    )
    if sort_by !== nothing
        sort!(tab_df, sort_by)
    end
    if compress === nothing
        compress = endswith(path, ".gz")
    end
    write_csv(path, tab_df; compress=compress)
    return nothing
end

function read_tabs(
    paths::Vector{String},
    key_col::Symbol,      # Name for the key column in output (e.g., :reads_per_umi)
    key_type::Type{K},    # Type for the key column in output (e.g., UInt64)
    value_col::Symbol,    # Name for the value column in output (e.g., :umis)
    value_type::Type{V},  # Type for the value column in output (e.g., Int64)
)::DataFrame where {K,V}
    tab_df = read_dfs(paths, Dict(key_col => key_type, value_col => value_type))
    print_start("Summing $value_col... ")
    tab_gdf = combine(groupby(tab_df, key_col), value_col => sum => value_col)
    println_done()
    return tab_gdf
end

function read_tabs_dict(
    paths::Vector{String},
    key_col::Symbol,      # Name for the key column in output (e.g., :reads_per_umi)
    ::Type{K},            # Type for the key column in output (e.g., UInt64)
    value_col::Symbol,    # Name for the value column in output (e.g., :umis)
    ::Type{V},            # Type for the value column in output (e.g., Int64)
)::Dict{K,V} where {K,V}
    tab_dict = Dict{K, V}()
    key_zero = zero(K)
    for path in paths
        print_start("Reading $(basename(path))... ")
        open(GzipDecompressorStream, path, "r") do io
            # validate that the header is correct
            header_line = readline(io)
            header_cols = split(header_line, ',')
            if header_cols[1] != string(key_col) || header_cols[2] != string(value_col)
                error("Unexpected header in $path: expected '$key_col,$value_col' but got '$header_line'")
            end
            keys = K[]
            values = V[]
            for line in eachline(io)
                cols = split(line, ',')
                push!(keys, parse(K, cols[1]))
                push!(values, parse(V, cols[2]))
            end
            # write_tab saves keys in the hash order of the Dict they came from. Inserting them
            # in that order into a smaller table piles them into adjacent slots and turns
            # probing quadratic (one 128 MB rpsb2 took 24 min vs ~20 s), so grow the table first
            sizehint!(tab_dict, length(tab_dict) + length(keys))
            for (key, value) in zip(keys, values)
                tab_dict[key] = get(tab_dict, key, key_zero) + value
            end
        end
        println_done()
    end
    return tab_dict
end
