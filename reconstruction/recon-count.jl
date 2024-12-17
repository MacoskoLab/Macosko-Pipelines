using CSV
using HDF5
using Plots
using Peaks: findminima
using ArgParse
using PDFmerger
using StatsBase
using StatsPlots
using Distributed
using KernelDensity
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
        "--regex", "-r"
        help = "Pattern to match FASTQ filenames"
        arg_type = String
        default = ".*"
        
        "--downsampling_level", "-p"
        help = "Level of downsampling"
        arg_type = Float64
        default = 1.0

        "--R1_barcodes", "-x"
        help = "Number of R1 barcodes (<1 means auto-pick)"
        arg_type = Int64
        default = 0

        "--R2_barcodes", "-y"
        help = "Number of R2 barcodes (<1 means auto-pick)"
        arg_type = Int64
        default = 0
    end

    return parse_args(ARGS, s)
end

# Load the command-line arguments
args = get_args()

const fastq_path = args["fastq_path"]
println("FASTQ path: "*fastq_path)
@assert isdir(fastq_path) "FASTQ path not found"
@assert !isempty(readdir(fastq_path)) "FASTQ path is empty"

const out_path = args["out_path"]
println("Output path: "*out_path)
Base.Filesystem.mkpath(out_path)
@assert isdir(out_path) "Output path could not be created"

const regex = Regex(args["regex"])
if regex != r".*"
    println("FASTQ regex: $regex")
end

const prob = args["downsampling_level"]
@assert 0 < prob <= 1 "Invalid downsampling level $prob"
if prob < 1
    println("Downsampling level: $prob")
end

const R1_barcodes = args["R1_barcodes"]
if R1_barcodes > 0
    println("R1 barcodes: $R1_barcodes")
else
    println("R1 barcodes: AUTO")
end

const R2_barcodes = args["R2_barcodes"]
if R2_barcodes > 0
    println("R2 barcodes: $R2_barcodes")
else
    println("R2 barcodes: AUTO")
end

println("Threads: $(Threads.nthreads())\n")

# Load the FASTQ paths
fastqs = readdir(fastq_path, join=true)
fastqs = filter(fastq -> endswith(fastq, ".fastq.gz"), fastqs)
fastqs = filter(fastq -> occursin(regex, fastq), fastqs)
@assert length(fastqs) >= 2 "ERROR: No FASTQ pairs found"
const R1s = filter(s -> occursin("_R1_", s), fastqs) ; println("R1s: ", basename.(R1s))
const R2s = filter(s -> occursin("_R2_", s), fastqs) ; println("R2s: ", basename.(R2s))
@assert length(R1s) > 0 && length(R2s) > 0 "ERROR: No FASTQ pairs found"
@assert length(R1s) == length(R2s) "ERROR: R1s and R2s are not all paired"
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s]
println("$(length(R1s)) pair(s) of FASTQs found\n")

# Create a plot showing the file paths
N = 30
p = plot(xlim=(0, 4), ylim=(0, N+1), framestyle=:none, size=(7*100, 8*100),
         legend=false, xticks=:none, yticks=:none)
annotate!(p, 0.1, N,   text("Input directory:", :left, 9))
annotate!(p, 0.1, N-1, text("$fastq_path", :left, 9))
annotate!(p, 0.1, N-2, text("Output directory:", :left, 9))
annotate!(p, 0.1, N-3, text("$out_path", :left, 9))
annotate!(p, 0.1, N-5, text("R1 FASTQs:", :left, 9))
annotate!(p, 2.1, N-5, text("R2 FASTQs:", :left, 9))
for (i, (R1, R2)) in enumerate(zip(R1s, R2s))
    i > N-6 && break
    annotate!(p, 0.1, N-5-i, text(basename(R1), :left, 9))
    annotate!(p, 2.1, N-5-i, text(basename(R2), :left, 9))
end
if length(R1s) > N-6
    annotate!(p, 2, 0, text("$(length(R1s)-N+6) FASTQ file(s) not shown", :center, 9))
end
savefig(p, joinpath(out_path, "filepaths.pdf"))

################################################################################

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
        @fastmath @inbounds u = [bases[(code >> n) & 3 + 1] for n in 0:2:16]
        return String(u)
    end
    
    @inline function encode_15(sb_1::SeqView, sb_2::SeqView)::UInt64
        @fastmath @inbounds b1 = dot(px8, (codeunits(sb_1) .>> 1) .& 3)
        @fastmath @inbounds b2 = dot(px7, (codeunits(sb_2) .>> 1) .& 3)
        return b1 + b2 * 4^8
    end
    @inline function decode_15(code::UInt64)::String
        @fastmath @inbounds u = [bases[(code >> n) & 3 + 1] for n in 0:2:28]
        return String(u)
    end
    
    @inline function encode_17(sb_1::SeqView, sb_2::SeqView)::UInt64
        @fastmath @inbounds b1 = dot(px9, (codeunits(sb_1) .>> 1) .& 3)
        @fastmath @inbounds b2 = dot(px8, (codeunits(sb_2) .>> 1) .& 3)
        return b1 + b2 * 4^9
    end
    @inline function decode_17(code::UInt64)::String
        @fastmath @inbounds u = [bases[(code >> n) & 3 + 1] for n in 0:2:32]
        return String(u)
    end

    # Determine the R1 bead type
    function learn_R1type(R1)
        iter = R1 |> open |> GzipDecompressorStream |> FASTQ.Reader
        s10 = 0 ; s17 = 0
        for (i, record) in enumerate(iter)
            i > 100000 ? break : nothing
            seq = FASTQ.sequence(record)
            length(seq) < 42 ? continue : nothing
            s10 += get_V10(seq)[3] == UP1
            length(seq) < 44 ? continue : nothing
            s17 += get_V17(seq)[3] == UP1
        end
        myid() == 1 && println("V10: ", s10, " V17: ", s17)
        return(s10 >= s17 ? "V10" : "V17")
    end
    
    R1_types = [learn_R1type(R1) for R1 in R1s]
    if all(x -> x == "V10", R1_types)
        const bead1_type = "V10"
        const R1_len = 42
        const get_R1 = get_V10
        const encode_sb1 = encode_15
        const decode_sb1 = decode_15
    elseif all(x -> x == "V17", R1_types)
        const bead1_type = "V17"
        const R1_len = 44
        const get_R1 = get_V17
        const encode_sb1 = encode_17
        const decode_sb1 = decode_17
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
            length(seq) < 34 ? continue : nothing
            s15 += get_V15(seq)[3] == UP2
            length(seq) < 36 ? continue : nothing
            s16 += get_V16(seq)[3] == UP2
        end
        myid() == 1 &&  println("V15: ", s15, " V16: ", s16)
        return(s15 >= s16 ? "V15" : "V16")
    end
    
    R2_types = [learn_R2type(R2) for R2 in R2s]
    if all(x -> x == "V15", R2_types)
        const bead2_type = "V15"
        const R2_len = 34
        const get_R2 = get_V15
        const encode_sb2 = encode_15
        const decode_sb2 = decode_15
    elseif all(x -> x == "V16", R2_types)
        const bead2_type = "V16"
        const R2_len = 36
        const get_R2 = get_V16
        const encode_sb2 = encode_17
        const decode_sb2 = decode_17
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
    const umi_homopolymer_whitelist = reduce(union, [listHDneighbors(c^9, i, bases) for c in bases for i in 0:2])
    const sbi_homopolymer_whitelist = Set(encode_str(str) for str in reduce(union, [listHDneighbors(c^15, i) for c in bases for i in 0:3]))
end

################################################################################

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
        end
    end
        
    return df, metadata
end

results = pmap(pair -> process_fastqs(pair...), zip(R1s, R2s))
rmprocs(workers())

df = vcat([r[1] for r in results]...)
metadata = reduce((x, y) -> mergewith(+, x, y), [r[2] for r in results])
results = nothing

println("...done") ; flush(stdout) ; GC.gc()

################################################################################

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

println("done") ; flush(stdout) ; GC.gc()
println("Total UMIs: $(nrow(df))") ; flush(stdout)
@assert nrow(df) > 0

################################################################################

# Save results of fastq parsing
# h5open(joinpath(out_path, "reads.h5"), "w") do file
#     file["sb1_2bit", compress=1] = df[!, :sb1_i]
#     file["umi1_2bit", compress=1] = df[!, :umi1_i]
#     file["sb2_2bit", compress=1] = df[!, :sb2_i]
#     file["umi2_2bit", compress=1] = df[!, :umi2_i]
#     file["reads", compress=1] = df[!, :reads]
# end

# Load previous fastq parsing results
# df = h5open(joinpath(out_path, "reads.h5"), "r") do file 
#     DataFrame(sb1_i = read(file["sb1_2bit"]),
#               umi1_i = read(file["umi1_2bit"]),
#               sb2_i = read(file["sb2_2bit"]),
#               umi2_i = read(file["umi2_2bit"]),
#               reads = read(file["reads"]))
# end
# metadata = Dict(String(row[1]) => parse(Int, row[2]) for row in CSV.Rows(joinpath(out_path, "metadata.csv"), header=false))


################################################################################

print("Computing barcode whitelist... ") ; flush(stdout)

# Helper method
function remove_intermediate(x, y)
    m = (y .!= vcat(y[2:end], NaN)) .| (y .!= vcat(NaN, y[1:end-1]))
    x = x[m] ; y = y[m]
    return(x, y)
end

# Count the number of times each barcode appears
const tab1 = countmap(df[!,:sb1_i])
const tab2 = countmap(df[!,:sb2_i])

# Automatic cutoff
# Use the elbow plot to determine which beads to use as our whitelist
#   The cutoff is auto-detected using the steepest part of the curve
#   To make finding it more consistent, set a reasonable min/max UMI cutoff
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

const uc1_auto = determine_umi_cutoff(tab1 |> values |> collect)
const uc2_auto = determine_umi_cutoff(tab2 |> values |> collect)
const bc1_auto = count(e -> e >= uc1_auto, tab1 |> values |> collect)
const bc2_auto = count(e -> e >= uc2_auto, tab2 |> values |> collect)

# Manual cutoff
# R1_barcodes: provided at command-line
# R2_barcodes: provided at command-line
function bc_to_uc(bc, table)
    kv = DataFrame(keys = collect(keys(table)), values = collect(values(table)))
    sort!(kv, :keys, rev=true)
    kv[!, :cumsum] = cumsum(kv[!, :values])
    uc = kv[!, :keys][findfirst(x -> x >= bc, kv[!, :cumsum])]
    return uc
end

const uc1_manual = R1_barcodes > 0 ? bc_to_uc(R1_barcodes, tab1 |> values |> countmap) : R1_barcodes
const uc2_manual = R2_barcodes > 0 ? bc_to_uc(R2_barcodes, tab2 |> values |> countmap) : R2_barcodes
const bc1_manual = R1_barcodes > 0 ? count(e -> e >= uc1_manual, tab1 |> values |> collect) : R1_barcodes
const bc2_manual = R2_barcodes > 0 ? count(e -> e >= uc2_manual, tab2 |> values |> collect) : R2_barcodes

# Use manual cutoff if > 0
if R1_barcodes > 0 
    const bc1 = bc1_manual
    const uc1 = uc1_manual
else
    const bc1 = bc1_auto
    const uc1 = uc1_auto
end
if R2_barcodes > 0 
    const bc2 = bc2_manual
    const uc2 = uc2_manual
else
    const bc2 = bc2_auto
    const uc2 = uc2_auto
end

# Plots
function umi_density_plot(table, uc_auto, uc_manual, R)
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
             title = "$R UMI Count Distribution", label = "Barcodes",
             titlefont=10, guidefont=8, legendfontsize=6, xlabelfontsize=8, ylabelfontsize=8)
    plot!(p, (10).^lx_s, (10).^ly_s, seriestype = :line, label="KDE")
    vline!(p, [uc_auto], linestyle = :dash, color = :red, label = "UMI cutoff (auto)")
    uc_manual > 0 && vline!(p, [uc_manual], linestyle = :dash, color = :purple, label = "UMI cutoff (manual)")
    xticks!(p, [10^i for i in 0:ceil(log10(maximum(x)))])
    yticks!(p, [10^i for i in 0:ceil(log10(maximum(y)))])
    return p
end

p1 = umi_density_plot(tab1 |> values |> countmap, uc1_auto, uc1_manual, "R1")
p3 = umi_density_plot(tab2 |> values |> countmap, uc2_auto, uc2_manual, "R2")

function elbow_plot(y, uc_auto, uc_manual, bc_auto, bc_manual, R)
    sort!(y, rev=true)
    x = 1:length(y)
    
    xp, yp = remove_intermediate(x, y)
    p = plot(xp, yp, seriestype = :line, xscale = :log10, yscale = :log10,
         xlabel = "$R Spatial Barcode Rank", ylabel = "UMI Count", 
         title = "$R Spatial Barcode Elbow Plot", label = "Barcodes",
         titlefont=10, guidefont=8, legendfontsize=6, xlabelfontsize=8, ylabelfontsize=8)
    hline!(p, [uc_auto], linestyle = :dash, color = :red, label = "UMI cutoff (auto)")
    uc_manual > 0 && hline!(p, [uc_manual], linestyle = :dash, color = :purple, label = "UMI cutoff (manual)")
    vline!(p, [bc_auto], linestyle = :dash, color = :green, label = "SB cutoff (auto)")
    bc_manual > 0 && vline!(p, [bc_manual], linestyle = :dash, color = :brown, label = "SB cutoff (manual)")
    xticks!(p, [10^i for i in 0:ceil(log10(maximum(xp)))])
    yticks!(p, [10^i for i in 0:ceil(log10(maximum(yp)))])
    return p
end

p2 = elbow_plot(tab1 |> values |> collect, uc1_auto, uc1_manual, bc1_auto, bc1_manual, "R1")
p4 = elbow_plot(tab2 |> values |> collect, uc2_auto, uc2_manual, bc2_auto, bc2_manual, "R2")

p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "elbows.pdf"))

metadata["R1_umicutoff_auto"] = uc1_auto
metadata["R2_umicutoff_auto"] = uc2_auto
metadata["R1_barcodes_auto"] = bc1_auto
metadata["R2_barcodes_auto"] = bc2_auto

metadata["R1_umicutoff_manual"] = uc1_manual
metadata["R2_umicutoff_manual"] = uc2_manual
metadata["R1_barcodes_manual"] = bc1_manual
metadata["R2_barcodes_manual"] = bc2_manual

metadata["R1_umicutoff"] = uc1
metadata["R2_umicutoff"] = uc2
metadata["R1_barcodes"] = bc1
metadata["R2_barcodes"] = bc2

println("done") ; flush(stdout) ; GC.gc()

################################################################################

print("Matching to barcode whitelist... ") ; flush(stdout)

function match_barcode(df, metadata)    
    wl1 = Set{UInt64}([k for (k, v) in tab1 if v >= uc1]) ; @assert length(wl1) == bc1
    wl2 = Set{UInt64}([k for (k, v) in tab2 if v >= uc2]) ; @assert length(wl2) == bc2

    m1 = [s1 in wl1 for s1 in df.sb1_i]
    m2 = [s2 in wl2 for s2 in df.sb2_i]

    metadata["R1_exact"] = sum(m1)
    metadata["R2_exact"] = sum(m2)

    df.keep = m1 .& m2
    filter!(:keep => identity, df)
    select!(df, Not(:keep))

    nothing
end
match_barcode(df, metadata) # this function modifies in-place
metadata["umis_exact"] = nrow(df)

println("done") ; flush(stdout) ; GC.gc()

################################################################################

print("Removing chimeras... ") ; flush(stdout)

function remove_chimeras(df, metadata)
    sort!(df, [:sb1_i, :umi1_i, :reads], rev = [false, false, true])
    before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb1_i,:umi1_i]]))
    after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb1_i,:umi1_i,:reads]]), false)
    df.chimeric1 = before_same .| after_same
    
    sort!(df, [:sb2_i, :umi2_i, :reads], rev = [false, false, true])
    before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb2_i,:umi2_i]]))
    after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb2_i,:umi2_i,:reads]]), false)
    df.chimeric2 = before_same .| after_same
    
    metadata["R1_chimeric"] = sum(df.chimeric1)
    metadata["R2_chimeric"] = sum(df.chimeric2)

    subset!(df, :chimeric1 => x -> .!x, :chimeric2 => x -> .!x)
    select!(df, Not([:chimeric1, :chimeric2]))
    nothing
end
remove_chimeras(df, metadata) # this function modifies in-place
metadata["umis_chimeric"] = metadata["umis_exact"] - nrow(df)

println("done") ; flush(stdout) ; GC.gc()

################################################################################

print("Counting UMIs... ") ; flush(stdout)

function hist_10cap(vec)
    tab = countmap(vec)
    plotdf = DataFrame(value = collect(keys(tab)), count = collect(values(tab)))
    sum10 = sum(filter(:value => v -> v >= 10, plotdf).count)
    filter!(:value => v -> v < 10, plotdf)
    push!(plotdf, (10, sum10))
    p = bar(plotdf.value, plotdf.count, legend = false,
            xticks = (1:10, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"]),
            titlefont = 10, guidefont = 8)
    return p
end
p1 = hist_10cap(df.reads)
xlabel!(p1, "Reads per UMI")
ylabel!(p1, "Number of filtered UMIs")
title!(p1, "Read depth")

function count_umis(df)
    select!(df, [:sb1_i, :sb2_i])
    sort!(df, [:sb1_i, :sb2_i])
    start = vcat(true, (df.sb1_i[2:end] .!= df.sb1_i[1:end-1]) .| (df.sb2_i[2:end] .!= df.sb2_i[1:end-1]))
    umis = vcat(diff(findall(start)), nrow(df)-findlast(start)+1)

    df.start = start ; filter!(:start => identity, df) ; select!(df, Not(:start))
    df.umi = umis
    nothing
end
count_umis(df) # this function modifies in-place

p2 = hist_10cap(df.umi)
xlabel!(p2, "UMIs per connection")
ylabel!(p2, "Number of connections")
title!(p2, "Connection distribution")

function sum_topn(v, n)
    return sum(sort(v, rev=true)[1:min(n, length(v))])
end
function plot_umi_distributions(df, col::Symbol)
    R = "R"*string(col)[3] # R1 or R2
    gdf = combine(groupby(df, col), :umi => sum => :umi,
                                    :umi => length => :connections,
                                    :umi => (v->sum_topn(v,5)) => :top5,
                                    :umi => (v->sum_topn(v,20)) => :top20,
                                    :umi => (v->sum_topn(v,50)) => :top50,
                                    :umi => (v->sum_topn(v,100)) => :top100)

    plotdf = vcat(DataFrame(x = 1, y = gdf.top5 ./ gdf.umi),
                  DataFrame(x = 2, y = gdf.top20 ./ gdf.umi),
                  DataFrame(x = 3, y = gdf.top50 ./ gdf.umi),
                  DataFrame(x = 4, y = gdf.top100 ./ gdf.umi))

    # SNR violins
    p1 = @df plotdf begin
        violin(:x, :y, line = 0, fill = (0.3, :blue), legend = false, titlefont = 10, guidefont = 8,
            xticks = ([1, 2, 3, 4], ["5", "20", "50", "100"]), yticks = [0.0, 0.25, 0.5, 0.75, 1.0],
        xlabel = "Number of top beads", ylabel = "%UMIs in top beads", title = "$R SNR")
        boxplot!(:x, :y, line = (1, :black), fill = (0.3, :grey), outliers = false, legend = false)
    end

    # log-umi vs log-connection density
    m = max(log10(maximum(gdf.umi)),log10(maximum(gdf.connections)))
    xmean = round(log10(mean(gdf.umi)), digits=2)
    xmed = round(log10(median(gdf.umi)), digits=2)
    ymean = round(log10(mean(gdf.connections)), digits=2)
    ymed = round(log10(median(gdf.connections)),digits=2)
    p2 = histogram2d(log10.(gdf.umi), log10.(gdf.connections),
            show_empty_bins=true, color=cgrad(:plasma, scale = :exp),
            xlabel="log10 UMIs (mean: $xmean, median: $xmed)",
            ylabel="log10 connections (mean: $ymean, median: $ymed)",
            title="$R UMI Distribution", titlefont = 10, guidefont = 8,
            xlims=(0, m), ylims=(0, m))
    plot!(p2, [0, m], [0, m], color=:black, linewidth=1, legend = false)
   
    return p1, p2
end
p3, p5 = plot_umi_distributions(df, :sb1_i)
p4, p6 = plot_umi_distributions(df, :sb2_i)

p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "SNR.pdf"))

p = plot(p5, p6, layout = (2, 1), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "histograms.pdf"))

println("done") ; flush(stdout) ; GC.gc()

################################################################################

print("Connection filter... ") ; flush(stdout)

function connection_filter(df::DataFrame, metadata::Dict, col::Symbol, z=+3)
    R = "R"*string(col)[3] # R1 or R2
    gdf = combine(groupby(df, col), :umi => sum => :umi,
                                    :umi => length => :connections,
                                    :umi => maximum => :maximum)

    logcon = log10.(gdf.connections)
    logmax = log10.(gdf.maximum)
    
    # assume values above mode follow a half-normal distribution
    kdens = kde(logcon)
    xmode = kdens.x[argmax(kdens.density)]
    sd = mean(filter(x -> x >= xmode, logcon) .- xmode) * sqrt(pi) / sqrt(2)
    z3 = xmode + sd * z

    # create the label
    a = round(10^z3, digits=2)
    b = sum(gdf.connections .> 10^z3)
    c = round(b / nrow(gdf) * 100, digits=2)
    label = "z = +$z: $a\nbeads: $b ($c%)"

    # connection plot
    p1 = barhist(logcon, bins=100, legend=false, line=0,
                 titlefont = 10, guidefont = 6, xlabelfontsize = 8, ylabelfontsize = 8,
                 title="$R Connections", xlabel="Connections (log10)", ylabel="Count")
    vline!(p1, [z3], color=:red, linestyle=:dash, linewidth=1)
    annotate!(p1, z3+Plots.xlims(p1)[2]*0.02, Plots.ylims(p1)[2]*0.95,
              text(label, :red, 6, :left))

    # max plot
    p2 = barhist(logmax, bins=100, legend=false, line=0,
                 titlefont = 10, guidefont = 6, xlabelfontsize = 8, ylabelfontsize = 8,
                 title="$R Max UMIs per connection", xlabel="Maximum UMI (log10)", ylabel="Count")
    
    gdf_remove = filter(:connections => c -> c > 10^z3, gdf)
    remove_set = Set(gdf_remove[:, col])

    metadata["$(R)_cxnfilter_z"] = z
    metadata["$(R)_cxnfilter_cutoff"] = ceil(10^z3)
    metadata["$(R)_cxnfilter_beads"] = length(remove_set)
    metadata["$(R)_cxnfilter"] = sum(gdf_remove.umi)
    
    return remove_set, p1, p2
end

R1_remove, p1, p3 = connection_filter(df, metadata, :sb1_i)
R2_remove, p2, p4 = connection_filter(df, metadata, :sb2_i)

p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "connection_filter.pdf"))

before = sum(df.umi)
subset!(df, :sb1_i => x -> .!in.(x, [R1_remove]), :sb2_i => y -> .!in.(y, [R2_remove]))
after = sum(df.umi)

metadata["umis_cxnfilter"] = before - after
metadata["umis_final"] = sum(df.umi)
metadata["connections_final"] = nrow(df)

println("done") ; flush(stdout) ; GC.gc()

################################################################################

# Compute more metadata
sequencing_saturation = round((1 - (metadata["umis_filtered"] / metadata["reads_filtered"]))*100, digits=1)
bead_ratio = round(metadata["R2_barcodes"]/metadata["R1_barcodes"], digits=2)
ubcf = metadata["umis_exact"] - metadata["umis_chimeric"]
metadata["R1_beadtype"] = parse(Int, join(filter(isdigit, bead1_type)))
metadata["R2_beadtype"] = parse(Int, join(filter(isdigit, bead2_type)))
metadata["downsampling_pct"] = round(Int, prob*100)

# Plot the metadata summary
m = metadata
function f(num)
    num = string(num)
    num = reverse(join([reverse(num)[i:min(i+2, end)] for i in 1:3:length(num)], ","))
    return(num)
end
function r(num1, num2)
    string(round(num1/num2*100, digits=2))*"%"
end
function d(num1, num2)
    f(num1)*" ("*r(num1, num2)*")"
end
data = [
("R1 bead type", "Total reads", "R2 bead type"),
(bead1_type, f(m["reads"]), bead2_type),
("R1 GG UP", prob<1 ? "Downsampling level" : "", "R2 GG UP"),
(d(m["R1_GG_UP"],m["reads"]), prob<1 ? "$prob" : "", d(m["R2_GG_UP"],m["reads"])),
("R1 no UP", "", "R2 no UP"),
(d(m["R1_no_UP"],m["reads"]), "", d(m["R2_no_UP"],m["reads"])),
("R1 LQ UMI" , "", "R2 LQ UMI"),
(d(m["R1_N_UMI"]+m["R1_homopolymer_UMI"],m["reads"]), "", d(m["R2_N_UMI"]+m["R2_homopolymer_UMI"],m["reads"])),
("R1 LQ SB", "", "R2 LQ SB"),    
(d(m["R1_N_SB"]+m["R1_homopolymer_SB"],m["reads"]), "", d(m["R2_N_SB"]+m["R2_homopolymer_SB"],m["reads"])),
("", "Filtered reads", ""),
("", d(m["reads_filtered"], m["reads"]), ""),
("", "Sequencing saturation", ""),
("", string(sequencing_saturation)*"%", ""),
("", "Filtered UMIs", ""),
(R1_barcodes > 0 ? "(manual)" : "", f(m["umis_filtered"]), R2_barcodes > 0 ? "(manual)" : ""),
("R1 UMI cutoff", "", "R2 UMI cutoff"),   
(f(m["R1_umicutoff"]), "", f(m["R2_umicutoff"])),
("R1 Barcodes", "R2:R1 ratio", "R2 Barcodes"),
(f(m["R1_barcodes"]), string(bead_ratio), f(m["R2_barcodes"])),
("R1 matched", "Matched UMIs", "R2 matched"),
(r(m["R1_exact"],m["umis_filtered"]), r(m["umis_exact"],m["umis_filtered"]), r(m["R2_exact"],m["umis_filtered"])),
("R1 chimeric", "Chimeric UMIs", "R2 chimeric"),
(r(m["R1_chimeric"],m["umis_exact"]), r(m["umis_chimeric"], m["umis_exact"]), r(m["R2_chimeric"],m["umis_exact"])),
("R1 high-cxn", "High-connection UMIs", "R2 high-cxn"),
(r(m["R1_cxnfilter"],ubcf), r(m["umis_cxnfilter"],ubcf), r(m["R2_cxnfilter"],ubcf)),
("", "Final UMIs", ""),
("", d(m["umis_final"],m["umis_filtered"]), ""),
]
p = plot(xlim=(0, 4), ylim=(0, 28+1), framestyle=:none, size=(7*100, 8*100),
         legend=false, xticks=:none, yticks=:none)
for (i, (str1, str2, str3)) in enumerate(data)
    annotate!(p, 1, 28 - i + 1, text(str1, :center, 12))
    annotate!(p, 2, 28 - i + 1, text(str2, :center, 12))
    annotate!(p, 3, 28 - i + 1, text(str3, :center, 12))
end
hline!(p, [14.5], linestyle = :solid, color = :black)
savefig(p, joinpath(out_path, "metadata.pdf"))

# Write the metadata to .csv
meta_df = DataFrame([Dict(:key => k, :value => v) for (k,v) in metadata])
sort!(meta_df, :key) ; meta_df = select(meta_df, :key, :value)
CSV.write(joinpath(out_path,"metadata.csv"), meta_df, writeheader=false)

################################################################################

print("Writing output... ") ; flush(stdout)

merge_pdfs([joinpath(out_path,"elbows.pdf"),
            joinpath(out_path,"metadata.pdf"),
            joinpath(out_path,"SNR.pdf"),
            joinpath(out_path,"connection_filter.pdf"),
            joinpath(out_path,"histograms.pdf"),
            joinpath(out_path,"filepaths.pdf")],
            joinpath(out_path,"QC.pdf"), cleanup=true)

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

# Save the matrix
using DelimitedFiles
open(GzipCompressorStream, joinpath(out_path, "sb1.txt.gz"), "w") do file
    writedlm(file, sb1_whitelist, "\n")
end
open(GzipCompressorStream, joinpath(out_path, "sb2.txt.gz"), "w") do file
    writedlm(file, sb2_whitelist, "\n")
end
rename!(df, Dict(:sb1_i => :sb1_index, :sb2_i => :sb2_index, :umi => :umi))
open(GzipCompressorStream, joinpath(out_path, "matrix.csv.gz"), "w") do file
    CSV.write(file, df, writeheader=true)
end

@assert all(f -> isfile(joinpath(out_path, f)), ["matrix.csv.gz", "sb1.txt.gz", "sb2.txt.gz", "QC.pdf", "metadata.csv"])

println("done") ; flush(stdout) ; GC.gc()
