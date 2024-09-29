using CSV
using Plots
using Peaks: findminima
using ArgParse
using PDFmerger
using StatsBase
using StatsPlots
using Distributed
using Distributions: pdf, Exponential

# R1 recognized bead types:
# JJJJJJJJ  TCTTCAGCGTTCCCGAGA JJJJJJJ  NNNNNNNVV (V10)
# JJJJJJJJJ TCTTCAGCGTTCCCGAGA JJJJJJJJ NNNNNNNNN (V17)
# R2 recognized bead types:
# JJJJJJJJJJJJJJJ   CTGTTTCCTG NNNNNNNNN          (V15)
# JJJJJJJJJJJJJJJJJ CTGTTTCCTG NNNNNNNNN          (V16)

# Read the command-line arguments
function get_args()
    s = ArgParseSettings()

    # Positional arguments
    @add_arg_table s begin
        "fastq_path"
        help = "Path to the directory of FASTQ files"
        arg_type = String
        required = true

        "out_path"
        help = "Output directory"
        arg_type = String
        required = true
    end
    
    # Optional arguments
    @add_arg_table s begin
        "--downsampling_level", "-p"
        help = "Level of downsampling"
        arg_type = Float64
        default = 1.0
    end

    return parse_args(ARGS, s)
end

# Load the command-line arguments
args = get_args()
const fastq_path = args["fastq_path"]
const out_path = args["out_path"]
const prob = args["downsampling_level"]

println("FASTQ path: "*fastq_path)
@assert isdir(fastq_path) "FASTQ path not found"
@assert !isempty(readdir(fastq_path)) "FASTQ path is empty"
println("Output path: "*out_path)
Base.Filesystem.mkpath(out_path)
@assert isdir(out_path) "Output path could not be created"
@assert 0 < prob <= 1 "Invalid downsampling level $prob"
if prob < 1
    println("Downsampling level: $prob")
end
println("Threads: $(Threads.nthreads())")

# Load the FASTQ paths
fastqs = readdir(fastq_path, join=true)
fastqs = filter(fastq -> endswith(fastq, ".fastq.gz"), fastqs)
@assert length(fastqs) > 1 "ERROR: No FASTQ pairs found"
const R1s = filter(s -> occursin("_R1_", s), fastqs) ; println("R1s: ", basename.(R1s))
const R2s = filter(s -> occursin("_R2_", s), fastqs) ; println("R2s: ", basename.(R2s))
@assert length(R1s) > 0 && length(R2s) > 0 "ERROR: No FASTQ pairs found"
@assert length(R1s) == length(R2s) "ERROR: R1s and R2s are not all paired"
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s]
println("$(length(R1s)) pair(s) of FASTQs found\n")

####################################################################################################

# Create a worker for each FASTQ pair
addprocs(length(R1s))

@everywhere begin
    using FASTX
    using CodecZlib
    using IterTools: product
    using DataFrames
    using StringViews
    using LinearAlgebra: dot
    using Combinatorics: combinations

    const R1s = $R1s
    const R2s = $R2s
    const prob = $prob

    # Read structure methods
    const SeqView = StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}}
    @inline function get_V10(seq::SeqView)
        @inbounds sb_1 = seq[1:8]
        @inbounds up = seq[9:26]
        @inbounds sb_2 = seq[27:33]
        @inbounds umi = seq[34:42]
        return sb_1, sb_2, up, umi
    end
    @inline function get_V17(seq::SeqView)
        @inbounds sb_1 = seq[1:9]
        @inbounds up = seq[10:27]
        @inbounds sb_2 = seq[28:35]
        @inbounds umi = seq[36:44]
        return sb_1, sb_2, up, umi
    end
    @inline function get_V15(seq::SeqView)
        @inbounds sb_1 = seq[1:8]
        @inbounds sb_2 = seq[9:15]
        @inbounds up = seq[16:25]
        @inbounds umi = seq[26:34]
        return sb_1, sb_2, up, umi
    end
    @inline function get_V16(seq::SeqView)
        @inbounds sb_1 = seq[1:9]
        @inbounds sb_2 = seq[10:17]
        @inbounds up = seq[18:27]
        @inbounds umi = seq[28:36]
        return sb_1, sb_2, up, umi
    end
    const UP1 = "TCTTCAGCGTTCCCGAGA"
    const UP2 = "CTGTTTCCTG"
    
    # String bit-encoding methods
    const bases = ['A','C','T','G'] # MUST NOT change this order
    const px7 = [convert(UInt32, 4^i) for i in 0:6]
    const px8 = [convert(UInt32, 4^i) for i in 0:7]
    const px9 = [convert(UInt32, 4^i) for i in 0:8]

    @inline function encode_str(str::String)::UInt64 # careful, encodes N as G
        return dot([4^i for i in 0:(length(str)-1)], (codeunits(str) .>> 1) .& 3)
    end
    
    @inline function encode_umi(umi::SeqView)::UInt32
        @fastmath @inbounds b = dot(px9, (codeunits(umi) .>> 1) .& 3)
        return b
    end
    @inline function decode_umi(code::UInt32)::String
        return String([bases[(code >> n) & 3 + 1] for n in 0:2:16])
    end
    
    @inline function encode_15(sb_1::SeqView, sb_2::SeqView)::UInt64
        @fastmath @inbounds b1 = dot(px8, (codeunits(sb_1) .>> 1) .& 3)
        @fastmath @inbounds b2 = dot(px7, (codeunits(sb_2) .>> 1) .& 3)
        return b1 + b2 * 4^8
    end
    @inline function decode_15(code::UInt64)::String
        return String([bases[(code >> n) & 3 + 1] for n in 0:2:28])
    end
    
    @inline function encode_17(sb_1::SeqView, sb_2::SeqView)::UInt64
        @fastmath @inbounds b1 = dot(px9, (codeunits(sb_1) .>> 1) .& 3)
        @fastmath @inbounds b2 = dot(px8, (codeunits(sb_2) .>> 1) .& 3)
        return b1 + b2 * 4^9
    end
    @inline function decode_17(code::UInt64)::String
        return String([bases[(code >> n) & 3 + 1] for n in 0:2:32])
    end

    # Determine the R1 bead type
    function learn_R1type(R1)
        iter = R1 |> open |> GzipDecompressorStream |> FASTQ.Reader
        s10 = 0 ; s17 = 0
        for (i, record) in enumerate(iter)
            i > 100000 ? break : nothing
            seq = FASTQ.sequence(record)
            length(seq) < 44 ? continue : nothing
            s10 += get_V10(seq)[3] == UP1
            s17 += get_V17(seq)[3] == UP1
        end
        myid() == 1 && println("V10: ", s10, " V17: ", s17)
        return(s10 >= s17 ? "V10" : "V17")
    end
    R1_types = [learn_R1type(R1) for R1 in R1s]
    if all(x -> x == "V10", R1_types)
        const bead1_type = "V10"
        const R1_len = 42
        get_R1 = get_V10
        encode_sb1 = encode_15
        decode_sb1 = decode_15
    elseif all(x -> x == "V17", R1_types)
        const bead1_type = "V17"
        const R1_len = 44
        get_R1 = get_V17
        encode_sb1 = encode_17
        decode_sb1 = decode_17
    else
        error("Error: The R1 bead type is not consistent ($R1_types)")
    end
    myid() == 1 && println("R1 bead type: $bead1_type")
    
    # Determine the R2 bead type
    function learn_R2type(R2)
        iter = R2 |> open |> GzipDecompressorStream |> FASTQ.Reader
        s15 = 0 ; s16 = 0
        for (i, record) in enumerate(iter)
            i > 100000 ? break : nothing
            seq = FASTQ.sequence(record)
            length(seq) < 36 ? continue : nothing
            s15 += get_V15(seq)[3] == UP2
            s16 += get_V16(seq)[3] == UP2
        end
        myid() == 1 &&  println("V15: ", s15, " V16: ", s16)
        return(s15 >= s16 ? "V15" : "V16")
    end
    R2_types = [learn_R2type(R2) for R2 in R2s]
    if all(x -> x == "V15", R2_types)
        const bead2_type = "V15"
        const R2_len = 34
        get_R2 = get_V15
        encode_sb2 = encode_15
        decode_sb2 = decode_15
    elseif all(x -> x == "V16", R2_types)
        const bead2_type = "V16"
        const R2_len = 36
        get_R2 = get_V16
        encode_sb2 = encode_17
        decode_sb2 = decode_17
    else
        error("Error: The R2 bead type is not consistent ($R2_types)")
    end
    myid() == 1 && println("R2 bead type: $bead2_type")
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
    
    const UP1_whitelist = reduce(union, [listHDneighbors(UP1, i) for i in 0:2])
    const UP2_whitelist = reduce(union, [listHDneighbors(UP2, i) for i in 0:1])
    const UP1_GG_whitelist = reduce(union, [listHDneighbors("G"^length(UP1), i) for i in 0:3])
    const UP2_GG_whitelist = reduce(union, [listHDneighbors("G"^length(UP2), i) for i in 0:2])
    const umi_homopolymer_whitelist = reduce(union, [listHDneighbors(c^9, i) for c in bases for i in 0:2])
    const sbi_homopolymer_whitelist = Set{UInt64}(encode_str(str) for str in reduce(union, [listHDneighbors(c^15, i) for c in bases for i in 0:3]))
end

####################################################################################################

println("\nReading FASTQs...") ; flush(stdout)

# Read the FASTQs
@everywhere function process_fastqs(R1, R2)
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
            break
        end
    end
        
    return df, metadata
end

@time results = pmap((pair) -> process_fastqs(pair...), zip(R1s, R2s))
df = vcat([view(r[1],:,:) for r in results]...)
metadata = reduce((x, y) -> mergewith(+, x, y), [r[2] for r in results])
results = nothing

rmprocs(workers())
println("...done") ; flush(stdout) ; GC.gc()

####################################################################################################
####################################################################################################

print("Counting reads... ") ; flush(stdout)

function count_reads(df, metadata)
    sort!(df, [:sb1_i, :sb2_i, :umi1_i, :umi2_i])
    @assert nrow(df) == metadata["reads_filtered"]
    start = vcat(true, reduce(.|, [df[2:end,c] .!= df[1:end-1,c] for c in names(df)]))
    df = df[start, :]
    df[!,:reads] = vcat(diff(findall(start)), metadata["reads_filtered"]-findlast(start)+1)
    @assert sum(df.reads) == metadata["reads_filtered"]
    metadata["umis_filtered"] = nrow(df)
    return df, metadata
end
df, metadata = count_reads(df, metadata)

# Save reads per umi distribution
function compute_rpu(df)
    rpu_dict = countmap(df[!,:reads])
    rpu_df = DataFrame(reads_per_umi = collect(keys(rpu_dict)), umis = collect(values(rpu_dict)))
    sort!(rpu_df, :reads_per_umi)
    return rpu_df
end
CSV.write(joinpath(out_path,"reads_per_umi.csv"), compute_rpu(df), writeheader=true)

println("done") ; flush(stdout) ; GC.gc()
println("Total UMIs: $(nrow(df))") ; flush(stdout)
@assert nrow(df) > 0

####################################################################################################

print("Computing barcode whitelist... ") ; flush(stdout)

function remove_intermediate(x, y)
    m = (y .!= vcat(y[2:end], NaN)) .| (y .!= vcat(NaN, y[1:end-1]))
    x = x[m] ; y = y[m]
    return(x, y)
end

# Use the elbow plot to determine which beads to use as our whitelist
#   The cutoff is auto-detected using the steepest part of the curve
#   To make finding it more consistent, we set a reasonable min/max UMI cutoff
#   uc (umi cutoff) is the steepest part of the curve between min_uc and max_uc
function determine_umi_cutoff(y)
    sort!(y, rev=true)
    x = 1:length(y)
    x, y = remove_intermediate(x, y)
    
    # find the steepest slope
    lx = log10.(x) ; ly = log10.(y)
    dydx = (ly[1:end-2] - ly[3:end]) ./ (lx[1:end-2] - lx[3:end])
    min_uc = 10 ; max_uc = 1000 ; m = log10(min_uc) .<= ly[2:end-1] .<= log10(max_uc)
    min_index = findall(m)[argmin(dydx[m])] + 1 + 2
    
    uc = round(Int64, 10^ly[min_index])
    return uc
end

const tab1 = countmap(df[!,:sb1_i])
const tab2 = countmap(df[!,:sb2_i])

const uc1 = determine_umi_cutoff(tab1 |> values |> collect)
const uc2 = determine_umi_cutoff(tab2 |> values |> collect)

const bc1 = count(e -> e >= uc1, tab1 |> values |> collect)
const bc2 = count(e -> e >= uc2, tab2 |> values |> collect)

function umi_density_plot(table, uc, R)
    x = collect(keys(table))
    y = collect(values(table))
    perm = sortperm(x)
    x = x[perm]
    y = y[perm]
    
    # Compute the KDE
    lx_s = 0:0.001:ceil(maximum(log10.(x)), digits=3)
    ly_s = []
    for lx_ in lx_s
        weights = [pdf(Exponential(0.05), abs(lx_ - lx)) for lx in log10.(x)]
        kde = sum(log10.(y) .* weights) / sum(weights)
        push!(ly_s, kde)
    end
    
    # Create a density plot
    p = plot(x, y, seriestype = :scatter, xscale = :log10, yscale = :log10, 
             xlabel = "Number of UMI", ylabel = "Frequency",
             markersize = 3, markerstrokewidth = 0.1,
             title = "$R UMI Count Distribution", titlefont = 10, guidefont = 8, label = "Barcodes")
    plot!(p, (10).^lx_s, (10).^ly_s, seriestype = :line, label="KDE")
    vline!(p, [uc], linestyle = :dash, color = :red, label = "UMI cutoff")
    xticks!(p, [10^i for i in 0:ceil(log10(maximum(x)))])
    yticks!(p, [10^i for i in 0:ceil(log10(maximum(y)))])
    return p
end

function elbow_plot(y, uc, bc, R)
    sort!(y, rev=true)
    x = 1:length(y)
    
    xp, yp = remove_intermediate(x, y)
    p = plot(xp, yp, seriestype = :line, xscale = :log10, yscale = :log10,
         xlabel = "$R Spatial Barcode Rank", ylabel = "UMI Count", 
         title = "$R Spatial Barcode Elbow Plot", titlefont = 10, guidefont = 8, label = "Barcodes")
    hline!(p, [uc], linestyle = :dash, color = :red, label = "UMI cutoff")
    vline!(p, [bc], linestyle = :dash, color = :green, label = "SB cutoff")
    xticks!(p, [10^i for i in 0:ceil(log10(maximum(xp)))])
    yticks!(p, [10^i for i in 0:ceil(log10(maximum(yp)))])
    return p
end

p1 = umi_density_plot(tab1 |> values |> countmap, uc1, "R1")
p3 = umi_density_plot(tab2 |> values |> countmap, uc2, "R2")

p2 = elbow_plot(tab1 |> values |> collect, uc1, bc1, "R1")
p4 = elbow_plot(tab2 |> values |> collect, uc2, bc2, "R2")

p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "elbows.pdf"))

metadata["R1_umicutoff"] = uc1
metadata["R2_umicutoff"] = uc2
metadata["R1_barcodes"] = bc1
metadata["R2_barcodes"] = bc2

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Matching to barcode whitelist... ") ; flush(stdout)

function match_barcode(df, tab1, tab2)
    matching_metadata = Dict("R1_exact"=>0, "R1_none"=>0, "R2_exact"=>0, "R2_none"=>0)
    
    wl1 = Set{UInt64}([k for (k, v) in tab1 if v >= uc1]) ; @assert length(wl1) == bc1
    wl2 = Set{UInt64}([k for (k, v) in tab2 if v >= uc2]) ; @assert length(wl2) == bc2

    m1 = [s1 in wl1 for s1 in df[!,:sb1_i]]
    m2 = [s2 in wl2 for s2 in df[!,:sb2_i]]

    matching_metadata["R1_exact"] = sum(m1)
    matching_metadata["R1_none"] = sum(.!m1)
    matching_metadata["R2_exact"] = sum(m2)
    matching_metadata["R2_none"] = sum(.!m2)

    df = df[m1 .& m2, :]

    return df, matching_metadata
end
df, matching_metadata = match_barcode(df, tab1, tab2)
metadata["umis_matched"] = nrow(df)

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Counting UMIs... ") ; flush(stdout)

function plot_reads_per_umi(vec)
    tab = countmap(vec)
    df = DataFrame(value = collect(keys(tab)), count = collect(values(tab)))
    sum10 = sum(df[df.value .> 10, :count])
    df = df[df.value .<= 10, :]
    idx = findfirst(df.value .== 10)
    if isnothing(idx)
        push!(df, (10, sum10))
    else
        df[idx, :count] += sum10
    end
    p = bar(df.value, df.count, legend = false,
            xlabel = "Reads per UMI", ylabel = "Number of filtered UMIs", title = "Read depth",
            xticks = (1:10, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"]),
            titlefont = 10, guidefont = 8)
    return p
end
p1 = plot_reads_per_umi(df[!,:reads])

function count_umis(df, metadata)
    # Remove chimeras
    sort!(df, [:sb1_i, :umi1_i, :reads], rev = [false, false, true])
    before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb1_i,:umi1_i]]))
    after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb1_i,:umi1_i,:reads]]), false)
    df.chimeric1 = before_same .| after_same
    
    sort!(df, [:sb2_i, :umi2_i, :reads], rev = [false, false, true])
    before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb2_i,:umi2_i]]))
    after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb2_i,:umi2_i,:reads]]), false)
    df.chimeric2 = before_same .| after_same
    
    metadata["umis_chimeric_R1"] = sum(df.chimeric1)
    metadata["umis_chimeric_R2"] = sum(df.chimeric2)

    subset!(df, :chimeric1 => x -> .!x, :chimeric2 => x -> .!x)
    select!(df, [:sb1_i, :sb2_i])

    # Count UMIs
    sort!(df, [:sb1_i, :sb2_i])
    bnd = vcat(true, (df.sb1_i[2:end] .!= df.sb1_i[1:end-1]) .| (df.sb2_i[2:end] .!= df.sb2_i[1:end-1]))
    umis = vcat(diff(findall(bnd)), nrow(df)-findlast(bnd)+1)
    df = df[bnd, :]
    df[!,:umi] = umis
    metadata["umis_final"] = sum(df.umi)
    metadata["connections_final"] = nrow(df)
        
    return df, metadata
end
df, metadata = count_umis(df, metadata)

function plot_umis_per_connection(vec)
    tab = countmap(vec)
    df = DataFrame(value = collect(keys(tab)), count = collect(values(tab)))
    sum10 = sum(df[df.value .> 10, :count])
    df = df[df.value .<= 10, :]
    idx = findfirst(df.value .== 10)
    if isnothing(idx)
        push!(df, (10, sum10))
    else
        df[idx, :count] += sum10
    end
    p = bar(df.value, df.count, legend = false,
            xlabel = "UMIs per connection", ylabel = "Number of connections", title = "Connection distribution",
            xticks = (1:10, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"]),
            titlefont = 10, guidefont = 8)
    return p
end
p2 = plot_umis_per_connection(df.umi)

# Factorize the barcode indexes
uniques1 = sort(collect(Set(df.sb1_i)))
uniques2 = sort(collect(Set(df.sb2_i)))
sb1_whitelist = [decode_sb1(sb1_i) for sb1_i in uniques1]
sb2_whitelist = [decode_sb2(sb2_i) for sb2_i in uniques2]
dict1 = Dict{UInt64, UInt64}(value => index for (index, value) in enumerate(uniques1))
dict2 = Dict{UInt64, UInt64}(value => index for (index, value) in enumerate(uniques2))
df.sb1_i = [dict1[k] for k in df.sb1_i]
df.sb2_i = [dict2[k] for k in df.sb2_i]
@assert sort(collect(Set(df.sb1_i))) == collect(1:length(Set(df.sb1_i)))
@assert sort(collect(Set(df.sb2_i))) == collect(1:length(Set(df.sb2_i)))

function sum_top5(v)
    return sum(sort(v, rev=true)[1:min(5, length(v))])
end
function sum_top20(v)
    return sum(sort(v, rev=true)[1:min(20, length(v))])
end
function sum_top50(v)
    return sum(sort(v, rev=true)[1:min(50, length(v))])
end
function sum_top100(v)
    return sum(sort(v, rev=true)[1:min(100, length(v))])
end
function plot_umi_distributions(df, col::Symbol)
    gdf = combine(groupby(df, col), :umi => sum => :umi,
                                    :umi => length => :count,
                                    :umi => maximum => :max,
                                    :umi => sum_top5 => :top5,
                                    :umi => sum_top20 => :top20,
                                    :umi => sum_top50 => :top50,
                                    :umi => sum_top100 => :top100)

    plotdf = vcat(DataFrame(x = 1, y = gdf.top5 ./ gdf.umi),
                  DataFrame(x = 2, y = gdf.top20 ./ gdf.umi),
                  DataFrame(x = 3, y = gdf.top50 ./ gdf.umi),
                  DataFrame(x = 4, y = gdf.top100 ./ gdf.umi))
    
    p1 = @df plotdf begin
        violin(:x, :y, line = 0, fill = (0.3, :blue), legend = false, titlefont = 10, guidefont = 8,
            xticks = ([1, 2, 3, 4], ["5", "20", "50", "100"]), yticks = [0.0, 0.25, 0.5, 0.75, 1.0],
        xlabel = "Number of top beads", ylabel = "%UMIs in top beads", title = "R$(string(col)[3]) SNR")
        boxplot!(:x, :y, line = (1, :black), fill = (0.3, :grey), outliers = false, legend = false)
    end

    m = max(log10(maximum(gdf.umi)),log10(maximum(gdf.count)))
    p2 = histogram2d(log10.(gdf.umi), log10.(gdf.count), show_empty_bins=true, color=cgrad(:plasma, scale = :exp),
            xlabel="log10 UMIs (mean: $(round(log10(mean(gdf.umi)),digits=2)), median: $(round(log10(median(gdf.umi)),digits=2)))",
            ylabel="log10 connections (mean: $(round(log10(mean(gdf.count)),digits=2)), median: $(round(log10(median(gdf.count)),digits=2)))",
            title="R$(string(col)[3]) UMI Distribution", titlefont = 10, guidefont = 8, xlims=(0, m), ylims=(0, m))
    plot!(p2, [0, m], [0, m], color=:black, linewidth=1, legend = false)

    select!(gdf, [col, :umi, :count, :max])
    rename!(gdf, :count => :connections)
    sort!(gdf, col)
    @assert gdf[!, col] == collect(1:nrow(gdf))
    @assert sum(gdf.umi) == sum(df.umi)
    @assert sum(gdf.connections) == nrow(df)
    
    return gdf, p1, p2
end

df1, p3, p5 = plot_umi_distributions(df, :sb1_i)
df2, p4, p6 = plot_umi_distributions(df, :sb2_i)

# Distrubtion of umi1 umi2, reads per umi, umi per connection
p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "SNR.pdf"))

p = plot(p5, p6, layout = (2, 1), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "histograms.pdf"))

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Writing output... ") ; flush(stdout)

# Write the metadata
function f(num)
    num = string(num)
    num = reverse(join([reverse(num)[i:min(i+2, end)] for i in 1:3:length(num)], ","))
    return(num)
end
function d(num1, num2)
    f(num1)*" ("*string(round(num1/num2*100, digits=2))*"%"*")"
end
m = metadata
mm = matching_metadata
data = [
("R1 bead type", "Total reads", "R2 bead type"),
(bead1_type, f(m["reads"]), bead2_type),
("R1 too short", prob<1 ? "Downsampling level" : "", "R2 too short"),
(d(m["R1_tooshort"],m["reads"]), prob<1 ? "$prob" : "", d(m["R2_tooshort"],m["reads"])),
("R1 GG UP", "", "R2 GG UP"),
(d(m["R1_GG_UP"],m["reads"]), "", d(m["R2_GG_UP"],m["reads"])),
("R1 no UP", "", "R2 no UP"),
(d(m["R1_no_UP"],m["reads"]), "", d(m["R2_no_UP"],m["reads"])),
("R1 LQ UMI" , "", "R2 LQ UMI"),
(d(m["R1_N_UMI"],m["reads"]), "", d(m["R2_N_UMI"],m["reads"])),
("R1 degen UMI", "", "R2 degen UMI"),
(d(m["R1_homopolymer_UMI"],m["reads"]), "", d(m["R2_homopolymer_UMI"],m["reads"])),
("R1 LQ SB", "", "R2 LQ SB"),
(d(m["R1_N_SB"],m["reads"]), "", d(m["R2_N_SB"],m["reads"])),
("R1 degen SB", "", "R2 degen SB"),
(d(m["R1_homopolymer_SB"],m["reads"]), "", d(m["R2_homopolymer_SB"],m["reads"])),
("", "Filtered reads", ""),
("", d(m["reads_filtered"], m["reads"]), ""),
("", "Sequencing saturation", ""),
("", string(round((1 - (m["umis_filtered"] / m["reads_filtered"]))*100, digits=1))*"%", ""),
("", "Filtered UMIs", ""),
("", f(m["umis_filtered"]), ""),
("R1 UMI cutoff: ", "", "R2 UMI cutoff"),
(f(m["R1_umicutoff"]), "", f(m["R2_umicutoff"])),
("R1 Barcodes", "R2:R1 ratio", "R2 Barcodes"),
(f(m["R1_barcodes"]), string(round(m["R2_barcodes"]/m["R1_barcodes"],digits=2)), f(m["R2_barcodes"])),
("R1 exact matches", "", "R2 exact matches"),
(d(mm["R1_exact"],m["umis_filtered"]), "", d(mm["R2_exact"],m["umis_filtered"])),
("", "Matched UMIs", ""),
("", d(m["umis_matched"], m["umis_filtered"]), ""),
("R1 chimeras", "", "R2 chimeras"),
(d(m["umis_chimeric_R1"],m["umis_matched"]), "", d(m["umis_chimeric_R2"],m["umis_matched"])),
("", "Final UMIs", ""),
("", d(m["umis_final"],m["umis_filtered"]), ""),
]
p = plot(xlim=(0, 4), ylim=(0, 34+1), framestyle=:none, size=(7*100, 8*100),
         legend=false, xticks=:none, yticks=:none)
for (i, (str1, str2, str3)) in enumerate(data)
    annotate!(p, 1, 34 - i + 1, text(str1, :center, 12))
    annotate!(p, 2, 34 - i + 1, text(str2, :center, 12))
    annotate!(p, 3, 34 - i + 1, text(str3, :center, 12))
end
hline!(p, [14.5], linestyle = :solid, color = :black)
savefig(p, joinpath(out_path, "metadata.pdf"))

merge_pdfs([joinpath(out_path,"elbows.pdf"),
            joinpath(out_path,"metadata.pdf"),
            joinpath(out_path,"SNR.pdf"),
            joinpath(out_path,"histograms.pdf")],
            joinpath(out_path,"QC.pdf"), cleanup=true)

# Save the metadata
metadata["R1_beadtype"] = parse(Int, join(filter(isdigit, bead1_type)))
metadata["R2_beadtype"] = parse(Int, join(filter(isdigit, bead2_type)))
metadata["downsampling_pct"] = round(Int, prob*100)
@assert isempty(intersect(keys(metadata), keys(matching_metadata)))
meta_df = DataFrame([Dict(:key => k, :value => v) for (k,v) in merge(metadata, matching_metadata)])
sort!(meta_df, :key) ; meta_df = select(meta_df, :key, :value)
CSV.write(joinpath(out_path,"metadata.csv"), meta_df, writeheader=false)

# Save the matrix
@assert length(df1.umi) == length(sb1_whitelist)
open(GzipCompressorStream, joinpath(out_path,"sb1.csv.gz"), "w") do file
    write(file, "sb1,umi,connections,max\n")
    for line in zip(sb1_whitelist, df1.umi, df1.connections, df1.max)
        write(file, join(line, ",") * "\n")
    end
end
@assert length(df2.umi) == length(sb2_whitelist)
open(GzipCompressorStream, joinpath(out_path,"sb2.csv.gz"), "w") do file
    write(file, "sb2,umi,connections,max\n")
    for line in zip(sb2_whitelist, df2.umi, df2.connections, df2.max)
        write(file, join(line, ",") * "\n")
    end
end
rename!(df, Dict(:sb1_i => :sb1_index, :sb2_i => :sb2_index, :umi => :umi))
open(GzipCompressorStream, joinpath(out_path,"matrix.csv.gz"), "w") do file
    CSV.write(file, df, writeheader=true)
end

@assert all(f -> isfile(joinpath(out_path, f)), ["matrix.csv.gz", "sb1.csv.gz", "sb2.csv.gz", "QC.pdf", "metadata.csv", "reads_per_umi.csv"])

println("done") ; flush(stdout) ; GC.gc()
