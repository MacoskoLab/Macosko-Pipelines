#NOTE: using Distributed and addprocs() should have already been called so that @everywhere and workers() work

# R1 recognized bead types:
# JJJJJJJJ  TCTTCAGCGTTCCCGAGA JJJJJJJ  NNNNNNNVV (V10)
# JJJJJJJJJ TCTTCAGCGTTCCCGAGA JJJJJJJJ NNNNNNNNN (V17)
# JJJJJJJJJ TCTTCAGCGT         JJJJJJJJ NNNNNNNNN (V19)
# R2 recognized bead types:
# JJJJJJJJJJJJJJJ   CTGTTTCCTG NNNNNNNNN          (V15)
# JJJJJJJJJJJJJJJJJ CTGTTTCCTG NNNNNNNNN          (V16)

@everywhere begin
    using HDF5
    using FASTX
    using CodecZlib
    using IterTools: product
    using DataFrames
    using LinearAlgebra: dot
    using Combinatorics: combinations
end

# Create fuzzy matching whitelists
@everywhere workers() begin
    function listHDneighbors(str, hd, charlist = ['A','C','G','T','N'])::Set{String}
        res = Set{String}()
        for inds in combinations(1:length(str), hd)
            chars = [str[i] for i in inds]
            pools = [setdiff(charlist, [char]) for char in chars]
            prods = product(pools...)
            for prod in prods
                s = str
                for (i, c) in zip(inds, prod)
                    s = s[1:i-1]*string(c)*s[i+1:end]
                end
                push!(res,s)
            end
        end
        return(res)
    end

    const UP1_whitelist = reduce(union, [listHDneighbors(up1, i) for i in 0:2 for up1 in [UP1, UP1[1:10]]])
    const UP2_whitelist = reduce(union, [listHDneighbors(UP2, i) for i in 0:1])
    const UP1_GG_whitelist = reduce(union, [listHDneighbors("G"^l, i) for i in 0:3 for l in [10,18]])
    const UP2_GG_whitelist = reduce(union, [listHDneighbors("G"^length(UP2), i) for i in 0:2])
    const umi_homopolymer_whitelist = reduce(union, [listHDneighbors(c^9, i, bases) for c in bases for i in 0:2])
    const sbi_homopolymer_whitelist = Set(encode_str(str) for str in reduce(union, [listHDneighbors(c^15, i) for c in bases for i in 0:3]))
end

# Read the FASTQs
@everywhere function process_fastqs(prob, bead1_type, bead2_type, R1, R2, tmp)
    bead1_info = bead1_type_to_info(bead1_type)
    bead2_info = bead2_type_to_info(bead2_type)
    R1_len = bead1_info.R_len
    get_R1 = bead1_info.get_R
    encode_sb1 = bead1_info.encode_sb
    decode_sb1 = bead1_info.decode_sb
    R2_len = bead2_info.R_len
    get_R2 = bead2_info.get_R
    encode_sb2 = bead2_info.encode_sb
    decode_sb2 = bead2_info.decode_sb
    it1 = R1 |> open |> GzipDecompressorStream |> FASTQ.Reader
    it2 = R2 |> open |> GzipDecompressorStream |> FASTQ.Reader

    df = DataFrame(sb1_i = UInt64[], umi1_i = UInt32[], sb2_i = UInt64[], umi2_i = UInt32[])
    metadata = Dict("reads"=>0, "reads_filtered"=>0,
                    "R1_tooshort"=>0, "R2_tooshort"=>0,
                    "R1_no_UP"=>0, "R2_no_UP"=>0, "R1_GG_UP"=>0, "R2_GG_UP"=>0,
                    "R1_N_UMI"=>0, "R2_N_UMI"=>0, "R1_homopolymer_UMI"=>0, "R2_homopolymer_UMI"=>0,
                    "R1_N_SB"=>0, "R2_N_SB"=>0, "R1_homopolymer_SB"=>0, "R2_homopolymer_SB"=>0)

    for record in zip(it1, it2)
        # Random dropout for downsampling
        prob < 1 && rand() > prob && continue

        metadata["reads"] += 1

        # Load the sequences
        seq1 = FASTQ.sequence(record[1])
        seq2 = FASTQ.sequence(record[2])

        # Validate the sequence length
        skip = false
        if length(seq1) < R1_len
            metadata["R1_tooshort"] += 1
            skip = true
        end
        if length(seq2) < R2_len
            metadata["R2_tooshort"] += 1
            skip = true
        end
        if skip
            continue
        end

        # Parse the read structure
        sb1_1, sb1_2, up1, umi1 = get_R1(seq1)
        sb2_1, sb2_2, up2, umi2 = get_R2(seq2)

        # Validate the UP
        skip = false
        if !in(up1, UP1_whitelist)
            metadata["R1_no_UP"] += 1
            skip = true
        end
        if !in(up2, UP2_whitelist)
            metadata["R2_no_UP"] += 1
            skip = true
        end
        if in(up1, UP1_GG_whitelist)
            metadata["R1_GG_UP"] += 1
            skip = true
        end
        if in(up2, UP2_GG_whitelist)
            metadata["R2_GG_UP"] += 1
            skip = true
        end
        if skip
            continue
        end

        # Validate the UMI
        skip = false
        if occursin('N', umi1)
            metadata["R1_N_UMI"] += 1
            skip = true
        end
        if occursin('N', umi2)
            metadata["R2_N_UMI"] += 1
            skip = true
        end
        if in(umi1, umi_homopolymer_whitelist)
            metadata["R1_homopolymer_UMI"] += 1
            skip = true
        end
        if in(umi2, umi_homopolymer_whitelist)
            metadata["R2_homopolymer_UMI"] += 1
            skip = true
        end
        if skip
            continue
        end

        # Check SB for N
        skip = false
        if occursin('N', sb1_1) || occursin('N', sb1_2)
            metadata["R1_N_SB"] += 1
            skip = true
        end
        if occursin('N', sb2_1) || occursin('N', sb2_2)
            metadata["R2_N_SB"] += 1
            skip = true
        end
        if skip
            continue
        end

        sb1_i = encode_sb1(sb1_1, sb1_2)
        sb2_i = encode_sb2(sb2_1, sb2_2)

        # Check SB for homopolymer
        skip = false
        if in(sb1_i & (4^15 - 1), sbi_homopolymer_whitelist)
            metadata["R1_homopolymer_SB"] += 1
            skip = true
        end
        if in(sb2_i & (4^15 - 1), sbi_homopolymer_whitelist)
            metadata["R2_homopolymer_SB"] += 1
            skip = true
        end
        if skip
            continue
        end

        # Update counts
        umi1_i = encode_umi(umi1)
        umi2_i = encode_umi(umi2)
        push!(df, (sb1_i, umi1_i, sb2_i, umi2_i))
        metadata["reads_filtered"] += 1

        if metadata["reads_filtered"] % 10_000_000 == 0
            println(metadata["reads_filtered"]) ; flush(stdout)
        end
    end

    # Write results to disk
    h5open(tmp, "w") do file
        file["sb1_i", compress=0] = df[!, :sb1_i]
        file["umi1_i", compress=0] = df[!, :umi1_i]
        file["sb2_i", compress=0] = df[!, :sb2_i]
        file["umi2_i", compress=0] = df[!, :umi2_i]
    end

    df = nothing ; GC.gc()
    return metadata
end

function process_results(out_path, prob, bead1_type, bead2_type, R1s, R2s)
    println("\nReading FASTQs...") ; flush(stdout)

    # Write the tmp files
    tmps = [joinpath(out_path, "tmp$(i).h5") for i in 1:length(R1s)]
    metadatas = pmap(pair -> process_fastqs(prob, bead1_type, bead2_type, pair...), zip(R1s, R2s, tmps))

    println("...done") ; flush(stdout) ; GC.gc()
    return metadatas
end

function combine_results(
    out_path::String,
    metadatas::Vector{Dict{String, Int64}},
)::Tuple{DataFrame, Dict{String, Int64}}
    print("Combining FASTQ results...") ; flush(stdout)

    tmps = [joinpath(out_path, "tmp$(i).h5") for i in 1:length(metadatas)]
    metadata = reduce((x, y) -> mergewith(+, x, y), metadatas)

    # Read the tmp files
    df = DataFrame(
        sb1_i  = Vector{UInt64}(undef, metadata["reads_filtered"]),
        umi1_i = Vector{UInt32}(undef, metadata["reads_filtered"]),
        sb2_i  = Vector{UInt64}(undef, metadata["reads_filtered"]),
        umi2_i = Vector{UInt32}(undef, metadata["reads_filtered"]),
    )
    pos = 1
    for i in 1:length(R1s)
        n = metadatas[i]["reads_filtered"]
        r = pos:(pos + n - 1)
        h5open(tmps[i], "r") do file
            df.sb1_i[r]  .= read(file["sb1_i"])
            df.umi1_i[r] .= read(file["umi1_i"])
            df.sb2_i[r]  .= read(file["sb2_i"])
            df.umi2_i[r] .= read(file["umi2_i"])
        end
        pos += n
    end
    @assert metadata["reads_filtered"] == pos - 1
    rm.(tmps)

    println("...done") ; flush(stdout) ; GC.gc()
    return df, metadata
end
