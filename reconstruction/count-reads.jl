using CSV
using StatsBase
using DataFrames

function count_reads(
    df::DataFrame,
    metadata::Dict{String,Int64},
    bead1_info::BeadInfo,
    bead2_info::BeadInfo,
    out_path::String,
)::Nothing
    decode_sb1 = bead1_info.decode_sb
    decode_sb2 = bead2_info.decode_sb

    print("Counting reads... ") ; flush(stdout)

    # Remove PCR duplicates
    function count_reads(df, metadata)
        sort!(df, [:sb1_i, :sb2_i, :umi1_i, :umi2_i])
        @assert nrow(df) == metadata["reads_filtered"]
        start = vcat(true, reduce(.|, [df[2:end,c] .!= df[1:end-1,c] for c in names(df)]))
        reads = vcat(diff(findall(start)), nrow(df)-findlast(start)+1)

        df.start = start ; filter!(:start => identity, df) ; select!(df, Not(:start))
        df.reads = reads
        @assert sum(df.reads) == metadata["reads_filtered"]
        metadata["umis_filtered"] = nrow(df)

        nothing
    end
    count_reads(df, metadata) # this function modifies in-place

    # # Save reads per umi distribution
    # function save_rpu(df, path)
    #     rpu_dict = countmap(df[!,:reads])
    #     rpu_df = DataFrame(reads_per_umi = collect(keys(rpu_dict)), umis = collect(values(rpu_dict)))
    #     sort!(rpu_df, :reads_per_umi)
    #     CSV.write(path, rpu_df, writeheader=true)
    # end
    # save_rpu(df, joinpath(out_path, "reads_per_umi.csv"))

    # # Save reads per bead, umis per bead
    # function save_rupb(df::DataFrame, col::Symbol, decode::Function, path::String)
    #     gdf = combine(groupby(df, col),
    #         :reads => sum => :reads,
    #         nrow => :umis
    #     )
    #     gdf[!,Symbol(String(col)[1:3])] = decode.(gdf[!,col])
    #     CSV.write(path, gdf[:, [4, 2, 3]], writeheader=true, compress=true)
    # end
    # save_rupb(df, :sb1_i, decode_sb1, joinpath(out_path, "readumi_per_sb1.csv.gz"))
    # save_rupb(df, :sb2_i, decode_sb2, joinpath(out_path, "readumi_per_sb2.csv.gz"))

    println("done") ; flush(stdout) ; GC.gc()
    println("Total UMIs: $(nrow(df))") ; flush(stdout)

    return nothing
end
