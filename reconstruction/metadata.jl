using CSV
using DataFrames

const all_metadata_keys = [
    "R1_GG_UP",
    "R1_N_SB",
    "R1_N_UMI",
    "R1_barcodes",
    "R1_barcodes_auto",
    "R1_barcodes_manual",
    "R1_beadtype",
    "R1_chimeric",
    "R1_cxnfilter",
    "R1_cxnfilter_beads",
    "R1_cxnfilter_cutoff",
    "R1_cxnfilter_z",
    "R1_exact",
    "R1_homopolymer_SB",
    "R1_homopolymer_UMI",
    "R1_no_UP",
    "R1_tooshort",
    "R1_umicutoff",
    "R1_umicutoff_auto",
    "R1_umicutoff_manual",
    "R2_GG_UP",
    "R2_N_SB",
    "R2_N_UMI",
    "R2_barcodes",
    "R2_barcodes_auto",
    "R2_barcodes_manual",
    "R2_beadtype",
    "R2_chimeric",
    "R2_cxnfilter",
    "R2_cxnfilter_beads",
    "R2_cxnfilter_cutoff",
    "R2_cxnfilter_z",
    "R2_exact",
    "R2_homopolymer_SB",
    "R2_homopolymer_UMI",
    "R2_no_UP",
    "R2_tooshort",
    "R2_umicutoff",
    "R2_umicutoff_auto",
    "R2_umicutoff_manual",
    "connections_final",
    "downsampling_pct",
    "reads",
    "reads_filtered",
    "umis_chimeric",
    "umis_cxnfilter",
    "umis_exact",
    "umis_filtered",
    "umis_final",
]

const unique_metadata_keys = [
    "R1_beadtype",
    "R2_beadtype",
    "downsampling_pct",
]

function write_metadata(out_path::String, metadata::Dict{String,Int64}, file_suffix::String = "")::Nothing
    if !all(k -> k in all_metadata_keys, keys(metadata))
        error("Metadata contains unrecognized keys: $(setdiff(keys(metadata), all_metadata_keys))")
    end
    metadata_df = DataFrame([Dict(:key => k, :value => v) for (k,v) in metadata])
    sort!(metadata_df, :key) ; metadata_df = select(metadata_df, :key, :value)
    CSV.write(joinpath(out_path,"metadata$file_suffix.csv"), metadata_df, writeheader=false)
    return nothing
end

function read_metadata_path(path::String)::Dict{String,Int64}
    metadata_df = CSV.read(path, DataFrame, header=false)
    rename!(metadata_df, Dict(:Column1 => :key, :Column2 => :value))
    metadata = Dict{String,Int64}(row.key => row.value for row in eachrow(metadata_df))
    if !all(k -> k in all_metadata_keys, keys(metadata))
        error("Metadata file $path contains unrecognized keys: $(setdiff(keys(metadata), all_metadata_keys))")
    end
    return metadata
end

function read_metadata(in_path::String, file_suffix::String = "")::Dict{String,Int64}
    print_start("Reading metadata... ")
    metadata = read_metadata_path(joinpath(in_path, "metadata$file_suffix.csv"))
    println_done()
    return metadata
end

function merge_metadata(metadata1::Dict{String,Int64}, metadata2::Dict{String,Int64})::Dict{String,Int64}
    merged = Dict{String,Int64}()
    for key in union(keys(metadata1), keys(metadata2))
        if key in unique_metadata_keys
            if haskey(metadata1, key) && haskey(metadata2, key)
                if metadata1[key] != metadata2[key]
                    error("Conflicting values for unique metadata key '$key': $(metadata1[key]) vs $(metadata2[key])")
                end
                merged[key] = metadata1[key]
            elseif haskey(metadata1, key)
                merged[key] = metadata1[key]
            else
                merged[key] = metadata2[key]
            end
        elseif key in all_metadata_keys
            merged[key] = get(metadata1, key, 0) + get(metadata2, key, 0)
        else
            if key in keys(metadata1) && key in keys(metadata2)
                error("Unrecognized metadata key '$key' found in both metadata dictionaries")
            elseif key in keys(metadata1)
                error("Unrecognized metadata key '$key' found in metadata1")
            else
                error("Unrecognized metadata key '$key' found in metadata2")
            end
        end
    end
    return merged
end

function find_metadatas(in_path::String)::Vector{String}
    metadatas = readdir(in_path, join=true)
    metadatas = filter(df -> startswith(basename(df), "metadata."), metadatas)
    metadatas = filter(df -> endswith(df, ".csv"), metadatas)
    println("metadatas: ", basename.(metadatas))
    @assert length(metadatas) > 0 "ERROR: No metadata CSVs found"
    println("$(length(metadatas)) metadata CSV(s) found\n")
    return metadatas
end

function read_metadatas(in_path::String)::Dict{String,Int64}
    metadatas = find_metadatas(in_path)
    print_start("Reading metadata... ")
    metadata = reduce(merge_metadata, [read_metadata_path(metadata_path) for metadata_path in metadatas])
    println_done()
    return metadata
end
