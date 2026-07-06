using CSV
using CodecZlib
using StatsBase
using DataFrames
using DataStructures

function count_reads(
    df::DataFrame,
    metadata::Dict{String,Int64},
    bead1_info::BeadInfo,
    bead2_info::BeadInfo,
    out_path::String,
    file_suffix::String = "",
)::Nothing
    decode_sb1 = bead1_info.decode_sb
    decode_sb2 = bead2_info.decode_sb

    print_start("Counting reads... ")

    # Remove PCR duplicates
    function count_reads(df, metadata)
        if nrow(df) == 0
            df.reads = Int64[]
            metadata["umis_filtered"] = 0
            return nothing
        end

        sort!(df, [:sb1_i, :sb2_i, :umi1_i, :umi2_i])
        @assert nrow(df) == metadata["reads_filtered"]
        start = vcat(true, reduce(.|, [df[2:end,c] .!= df[1:end-1,c] for c in names(df)]))
        reads = vcat(diff(findall(start)), nrow(df)-findlast(start)+1)

        df.start = start ; filter!(:start => identity, df) ; select!(df, Not(:start))
        df.reads = reads
        @assert sum(df.reads) == metadata["reads_filtered"]
        metadata["umis_filtered"] = nrow(df)

        return nothing
    end
    count_reads(df, metadata) # this function modifies in-place

    # # Save reads per umi distribution
    # function save_rpu(df, path)
    #     write_tab(df, path, :reads, :reads_per_umi, :umis; sort_by=:reads_per_umi)
    # end
    # save_rpu(df, joinpath(out_path, "reads_per_umi$file_suffix.csv"))

    # # Save reads per bead, umis per bead
    # function save_rupb(df::DataFrame, col::Symbol, decode::Function, path::String)
    #     gdf = combine(groupby(df, col),
    #         :reads => sum => :reads,
    #         nrow => :umis
    #     )
    #     gdf[!,Symbol(String(col)[1:3])] = decode.(gdf[!,col])
    #     CSV.write(path, gdf[:, [4, 2, 3]], writeheader=true, compress=true)
    # end
    # save_rupb(df, :sb1_i, decode_sb1, joinpath(out_path, "readumi_per_sb1$file_suffix.csv.gz"))
    # save_rupb(df, :sb2_i, decode_sb2, joinpath(out_path, "readumi_per_sb2$file_suffix.csv.gz"))

    println_done()
    println("Total UMIs: $(nrow(df))") ; flush(stdout)

    return nothing
end

function find_rpus(in_path::String, file_suffix::String)::Vector{String}
    rpus = readdir(in_path, join=true)
    rpus = filter(rpu -> startswith(basename(rpu), "reads_per_umi."), rpus)
    rpus = filter(rpu -> endswith(rpu, ".csv"), rpus)
    if length(file_suffix) > 0
        rpus = filter(rpu -> contains(rpu, file_suffix), rpus)
    end
    println("rpus: ", basename.(rpus))
    @assert length(rpus) > 0 "ERROR: No reads per umi CSVs found"
    println("$(length(rpus)) reads per umi CSV(s) found\n")
    return rpus
end

function merge_rpus(in_path::String, out_path::String, file_suffix::String)::Nothing
    rpus = find_rpus(in_path, file_suffix)
    df = read_tabs(rpus, :reads_per_umi, Int64, :umis, Int64)
    sort_df(df, :reads_per_umi)
    write_df(df, joinpath(out_path, "reads_per_umi.csv"))
    return nothing
end

function find_rupbs(in_path::String, col::Symbol, file_suffix::String)::Vector{String}
    readumi_per_sbs = readdir(in_path, join=true)
    readumi_per_sbs = filter(
        readumi_per_sb -> startswith(basename(readumi_per_sb), "readumi_per_$col."),
        readumi_per_sbs,
    )
    readumi_per_sbs = filter(readumi_per_sb -> endswith(readumi_per_sb, ".csv.gz"), readumi_per_sbs)
    if length(file_suffix) > 0
        readumi_per_sbs = filter(readumi_per_sb -> contains(readumi_per_sb, file_suffix), readumi_per_sbs)
    end
    println("readumi_per_sb1s: ", basename.(readumi_per_sbs))
    @assert length(readumi_per_sbs) > 0 "ERROR: No readumi per $col CSVs found"
    println("$(length(readumi_per_sbs)) of readumi per $col found\n")
    return readumi_per_sbs
end

function read_rupbs_dict(paths::Vector{String}, col::Symbol)::SortedDict{String,Tuple{UInt64,UInt64}}
    rupbs_dict = SortedDict{String, Tuple{UInt64,UInt64}}()
    for path in paths
        print_start("Reading $(basename(path))... ")
        open(GzipDecompressorStream, path, "r") do io
            header_line = readline(io)
            if header_line != "$col,reads,umis"
                error("Unexpected header in $path: expected '$col,reads,umis' but got '$header_line'")
            end
            for line in eachline(io)
                parts = split(line, ',')
                key = parts[1]
                reads = parse(UInt64, parts[2])
                umis = parse(UInt64, parts[3])
                pair = get(rupbs_dict, key, (0, 0))
                rupbs_dict[key] = (pair[1] + reads, pair[2] + umis)
            end
        end
        println_done()
    end
    return rupbs_dict
end

function write_rupbs_dict(
    rupbs_dict::SortedDict{String,Tuple{UInt64,UInt64}},
    out_path::String,
    col::Symbol,
    file_suffix::String,
)::Nothing
    print_start("Writing readumi per $col... ")
    open(GzipCompressorStream, joinpath(out_path, "readumi_per_$col$file_suffix.csv.gz"), "w") do io
        println(io, "$col,reads,umis")
        for (key, (reads, umis)) in rupbs_dict
            println(io, "$key,$reads,$umis")
        end
    end
    println_done()
    return nothing
end

function merge_rupbs1(in_path::String, out_path::String, file_suffix::String = "")::Nothing
    readumi_per_sbs = find_rupbs(in_path, :sb1, file_suffix)
    rupbs_dict = read_rupbs_dict(readumi_per_sbs, :sb1)
    write_rupbs_dict(rupbs_dict, out_path, :sb1, file_suffix)
    return nothing
end

function merge_rupbs2(in_path::String, out_path::String, file_suffix::String = "")::Nothing
    readumi_per_sbs = find_rupbs(in_path, :sb2, file_suffix)
    rupbs_dict = read_rupbs_dict(readumi_per_sbs, :sb2)
    write_rupbs_dict(rupbs_dict, out_path, :sb2, file_suffix)
    return nothing
end

function find_rpsbs(in_path::String)::Tuple{Vector{String}, Vector{String}}
    rpsbs = readdir(in_path, join=true)
    rpsbs = filter(rpsb -> startswith(basename(rpsb), "rpsb1.") || startswith(basename(rpsb), "rpsb2."), rpsbs)
    rpsbs = filter(rpsb -> endswith(rpsb, ".csv.gz"), rpsbs)
    @assert length(rpsbs) >= 2 "ERROR: No reads per sb pairs found"
    rpsb1s = filter(s -> startswith(basename(s), "rpsb1."), rpsbs) ; println("rpsb1s: ", basename.(rpsb1s)) # const
    rpsb2s = filter(s -> startswith(basename(s), "rpsb2."), rpsbs) ; println("rpsb2s: ", basename.(rpsb2s)) # const
    @assert length(rpsb1s) > 0 && length(rpsb2s) > 0 "ERROR: No reads per sb pairs found"
    @assert length(rpsb1s) == length(rpsb2s) "ERROR: rpsb1s and rpsb2s are not all paired"
    println("$(length(rpsb1s)) pair(s) of reads per sbs found\n")
    return rpsb1s, rpsb2s
end

function read_rpsbs(in_path::String)::Tuple{Dict{UInt64,Int64}, Dict{UInt64,Int64}}
    rpsb1s, rpsb2s = find_rpsbs(in_path)
    rpsb1 = read_tabs_dict(rpsb1s, :sb1_i, UInt64, :reads, Int64)
    rpsb2 = read_tabs_dict(rpsb2s, :sb2_i, UInt64, :reads, Int64)
    return rpsb1, rpsb2
end

function write_rpsbs(df::DataFrame, out_path::String, file_suffix::String)::Nothing
    print_start("Writing reads per sb1... ")
    write_tab(df, joinpath(out_path, "rpsb1$file_suffix.csv.gz"), :sb1_i, :sb1_i, :reads)
    println_done()
    print_start("Writing reads per sb2... ")
    write_tab(df, joinpath(out_path, "rpsb2$file_suffix.csv.gz"), :sb2_i, :sb2_i, :reads)
    println_done()
    return nothing
end

function find_reads_dfs(in_path::String, file_suffix::String = "")::Vector{String}
    dfs = readdir(in_path, join=true)
    dfs = filter(df -> startswith(basename(df), "reads_df."), dfs)
    dfs = filter(df -> endswith(df, ".csv.gz"), dfs)
    if length(file_suffix) > 0
        dfs = filter(df -> contains(df, file_suffix), dfs)
    end
    println("reads_dfs: ", basename.(dfs))
    @assert length(dfs) > 0 "ERROR: No reads data frames found"
    println("$(length(dfs)) reads data frame(s) found\n")
    return dfs
end

function read_reads_dfs(in_path::String, file_suffix::String = "")::DataFrame
    types = Dict(:sb1_i => UInt64, :sb2_i => UInt64, :umi1_i => UInt64, :umi2_i => UInt64, :reads => Int64)
    paths = find_reads_dfs(in_path, file_suffix)
    return read_dfs(paths, types)
end

function write_reads_df(df::DataFrame, out_path::String, file_suffix::String)::Nothing
    write_df(df, joinpath(out_path, "reads_df$file_suffix.csv.gz"))
    return nothing
end

function read_reads_df(in_path::String, file_suffix::String)::DataFrame
    types = Dict(:sb1_i => UInt64, :sb2_i => UInt64, :umi1_i => UInt64, :umi2_i => UInt64, :reads => Int64)
    return read_df(joinpath(in_path, "reads_df$file_suffix.csv.gz"), types)
end
