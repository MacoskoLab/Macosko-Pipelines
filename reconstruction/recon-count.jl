using CSV
using FASTX
using Plots
using Loess
using Peaks: findminima
using CodecZlib
using PDFmerger
using IterTools: product
using StatsBase: countmap
using DataFrames
using StringViews
using LinearAlgebra: dot
using Combinatorics: combinations

# Load the command-line arguments
if length(ARGS) != 2
    error("Usage: julia recon-count.jl fastq_path out_path")
end
fastq_path = ARGS[1]
println("FASTQ path: "*fastq_path)
@assert isdir(fastq_path) "FASTQ path not found"
@assert !isempty(readdir(fastq_path)) "FASTQ path is empty"
out_path = ARGS[2]
println("Output path: "*out_path)
Base.Filesystem.mkpath(out_path)
@assert isdir(out_path) "Output path not created"

# Load the FASTQ paths
fastqs = readdir(fastq_path, join=true)
R1s = filter(s -> occursin("_R1_", s), fastqs) ; println("R1s: ", basename.(R1s))
R2s = filter(s -> occursin("_R2_", s), fastqs) ; println("R2s: ", basename.(R2s))
@assert length(R1s) == length(R2s) > 0
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s]

# Validate the FASTQ sequence lengths
function fastq_seq_len(path)
    return(path |> open |> GzipDecompressorStream |> FASTQ.Reader |> first |> FASTQ.sequence |> length)
end
@assert length(unique([fastq_seq_len(R1) for R1 in R1s])) == 1 "WARNING: R1s have different FASTQ sequence lengths - proceed only if you are sure they have the same read structure"
@assert length(unique([fastq_seq_len(R2) for R2 in R2s])) == 1 "WARNING: R2s have different FASTQ sequence lengths - proceed only if you are sure they have the same read structure"

# UP matching methods
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
const UP1 = "TCTTCAGCGTTCCCGAGA"
const UP2 = "CTGTTTCCTG"
const UP1_whitelist = Set{String31}(reduce(union, [listHDneighbors(UP1, i) for i in 0:3]))
const UP2_whitelist = Set{String15}(reduce(union, [listHDneighbors(UP2, i) for i in 0:2]))
const umi_homopolymer_whitelist = Set{String15}(reduce(union, [listHDneighbors(str, i) for str in [c^9 for c in ["A","C","G","T"]] for i in 0:2]))

# UMI compressing methods
const px = [convert(UInt32, 4^i) for i in 0:(9-1)]
function UMItoindex(UMI::StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}})::UInt32
    return(dot(px, (codeunits(UMI).>>1).&3))
end
const bases = ['A','C','T','G'] # MUST NOT change this order
function indextoUMI(i::UInt32)::String15
    return(String15(String([bases[(i>>n)&3+1] for n in 0:2:16])))
end

println("R1 bead type: V8/V10")
println("R2 bead type: V15")

####################################################################################################

print("Reading FASTQs... ") ; flush(stdout)

# Read the FASTQs
function process_fastqs(R1s, R2s)
    sb1_dictionary = Dict{String15, UInt32}() # sb1 -> sb1_i
    sb2_dictionary = Dict{String15, UInt32}() # sb2 -> sb2_i
    mat = Dict{Tuple{UInt32, UInt32, UInt32, UInt32}, UInt32}() # (sb1_i, umi1_i, sb2_i, umi2_i) -> reads
    metadata = Dict("reads"=>0, "R1_tooshort"=>0, "R2_tooshort"=>0, "R1_N_UMI"=>0, "R2_N_UMI"=>0, "R1_homopolymer_UMI"=>0, "R2_homopolymer_UMI"=>0, "R1_no_UP"=>0, "R2_no_UP"=>0, "reads_filtered"=>0)

    for fastqpair in zip(R1s, R2s)
        it1 = fastqpair[1] |> open |> GzipDecompressorStream |> FASTQ.Reader
        it2 = fastqpair[2] |> open |> GzipDecompressorStream |> FASTQ.Reader
        for record in zip(it1, it2)
            metadata["reads"] += 1

            skip = false
            if length(FASTQ.sequence(record[1])) < 42
                metadata["R1_tooshort"] += 1
                skip = true
            end
            if length(FASTQ.sequence(record[2])) < 34
                metadata["R2_tooshort"] += 1
                skip = true
            end
            if skip
                continue
            end

            sb1_1 = FASTQ.sequence(record[1], 1:8)
            up1 = FASTQ.sequence(record[1], 9:26)
            sb1_2 = FASTQ.sequence(record[1], 27:33)
            umi1 = FASTQ.sequence(record[1], 34:42)

            sb2_1 = FASTQ.sequence(record[2], 1:8)
            sb2_2 = FASTQ.sequence(record[2], 9:15)
            up2 = FASTQ.sequence(record[2], 16:25)
            umi2 = FASTQ.sequence(record[2], 26:34)

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
            if !in(up1, UP1_whitelist)
                metadata["R1_no_UP"] += 1
                skip = true
            end
            if !in(up2, UP2_whitelist)
                metadata["R2_no_UP"] += 1
                skip = true
            end
            if skip
                continue
            end

            # update counts
            sb1_i = get!(sb1_dictionary, sb1_1*sb1_2, length(sb1_dictionary) + 1)
            sb2_i = get!(sb2_dictionary, sb2_1*sb2_2, length(sb2_dictionary) + 1)
            umi1_i = UMItoindex(umi1)
            umi2_i = UMItoindex(umi2)
            key = (sb1_i, umi1_i, sb2_i, umi2_i)
            mat[key] = get(mat, key, 0) + 1
            metadata["reads_filtered"] += 1
        end
    end

    sb1_whitelist = DataFrame(sb1 = collect(String15, keys(sb1_dictionary)), sb1_i = collect(UInt32, values(sb1_dictionary)))
    sb2_whitelist = DataFrame(sb2 = collect(String15, keys(sb2_dictionary)), sb2_i = collect(UInt32, values(sb2_dictionary)))
    sort!(sb1_whitelist, :sb1_i)
    sort!(sb2_whitelist, :sb2_i)
    @assert sb1_whitelist.sb1_i == 1:size(sb1_whitelist, 1)
    @assert sb2_whitelist.sb2_i == 1:size(sb2_whitelist, 1)
    
    metadata["umis_filtered"] = length(mat)
    return(mat, sb1_whitelist.sb1, sb2_whitelist.sb2, metadata)
end
mat, sb1_whitelist, sb2_whitelist, metadata = process_fastqs(R1s, R2s)

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Computing barcode whitelist... ") ; flush(stdout)

# Create plots and determine cutoff
function umi_density_plot(table, R)
    x = collect(keys(table))
    y = collect(values(table))
    perm = sortperm(x)
    x = x[perm]
    y = y[perm]
    
    model = loess(log10.(x), log10.(y), span=0.075)
    ys = (10).^predict(model, log10.(x))
    uc = x[findminima(ys).indices]
    uc = filter(x -> x >= 10, uc) |> sort |> first
    
    p = plot(x, y, seriestype = :scatter, xscale = :log10, yscale = :log10, 
             xlabel = "Number of UMI", ylabel = "Frequency",
             markersize = 3, markerstrokewidth = 0.1,
             title = "$R UMI Count Distribution", titlefont = 10, guidefont = 8, label = "Barcodes")
    plot!(p, x, ys, seriestype = :line, label="LOESS")
    vline!(p, [uc], linestyle = :dash, color = :red, label = "UMI cutoff")
    xticks!(p, [10^i for i in 0:ceil(log10(maximum(x)))])
    yticks!(p, [10^i for i in 0:ceil(log10(maximum(y)))])
    return(p, uc)
end
function remove_int(x, y)
    m = (y .!= vcat(y[2:end], NaN)) .| (y .!= vcat(NaN, y[1:end-1]))
    x = x[m]
    y = y[m]
    return(x, y)
end
function elbow_plot(y, uc, R)
    sort!(y, rev=true)
    x = 1:length(y)
    bc = count(e -> e >= uc, y)
    
    xp, yp = remove_int(x, y)
    p = plot(xp, yp, seriestype = :line, xscale = :log10, yscale = :log10,
         xlabel = "$R Spatial Barcode Rank", ylabel = "UMI Count", 
         title = "$R Spatial Barcode Elbow Plot", titlefont = 10, guidefont = 8, label = "Barcodes")
    hline!(p, [uc], linestyle = :dash, color = :red, label = "UMI cutoff")
    vline!(p, [bc], linestyle = :dash, color = :red, label = "SB cutoff")
    xticks!(p, [10^i for i in 0:ceil(log10(maximum(xp)))])
    yticks!(p, [10^i for i in 0:ceil(log10(maximum(yp)))])
    return(p, bc)
end

tab1 = countmap([key[1] for key in keys(mat)])
tab2 = countmap([key[3] for key in keys(mat)])

p1, uc1 = umi_density_plot(tab1 |> values |> countmap, "R1")
p3, uc2 = umi_density_plot(tab2 |> values |> countmap, "R2")

p2, bc1 = elbow_plot(tab1 |> values |> collect, uc1, "R1")
p4, bc2 = elbow_plot(tab2 |> values |> collect, uc2, "R2")

p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
savefig(p, joinpath(out_path, "elbows.pdf"))

metadata["R1_umicutoff"] = uc1
metadata["R2_umicutoff"] = uc2
metadata["R1_barcodes"] = bc1
metadata["R2_barcodes"] = bc2

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Matching to barcode whitelist... ") ; flush(stdout)

wl1 = Set{UInt32}([k for (k, v) in tab1 if v >= uc1]) ; @assert length(wl1) == bc1
wl2 = Set{UInt32}([k for (k, v) in tab2 if v >= uc2]) ; @assert length(wl2) == bc2

function match_barcode(mat, wl1, wl2)
    matching_metadata = Dict("R1_exact"=>0, "R1_none"=>0,
                             "R2_exact"=>0, "R2_none"=>0)
    
    for key in keys(mat)
        if key[1] in wl1
            matching_metadata["R1_exact"] += 1
        else
            matching_metadata["R1_none"] += 1
            delete!(mat, key)
        end
        if key[3] in wl2
            matching_metadata["R2_exact"] += 1
        else
            matching_metadata["R2_none"] += 1
            delete!(mat, key)
        end
    end
    return mat, matching_metadata
end
mat, matching_metadata = match_barcode(mat, wl1, wl2)
metadata["umis_matched"] = length(mat)

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Counting UMIs... ") ; flush(stdout)

function count_umis(mat)
    umi_dict = Dict{Tuple{UInt32, UInt32}, UInt32}()
    for key in keys(mat)
        pair = (key[1], key[3])
        pop!(mat, key)
        umi_dict[pair] = get(umi_dict, pair, 0) + 1
    end
    df = DataFrame(sb1_i = UInt32[], sb2_i = UInt32[], umi = UInt32[])
    for key in keys(umi_dict)
        value = pop!(umi_dict, key)
        push!(df, (key[1], key[2], value))
    end
    return(df)
end
df = count_umis(mat)

println("done") ; flush(stdout) ; GC.gc()

####################################################################################################

print("Writing Output... ") ; flush(stdout)

# Write the metadata
s = string
function f(num)
    num = string(num)
    num = reverse(join([reverse(num)[i:min(i+2, end)] for i in 1:3:length(num)], ","))
    return(num)
end
function d(num1, num2)
    f(num1)*" ("*s(round(num1/num2*100, digits=2))*"%"*")"
end
m = metadata
mm = matching_metadata
data = [
("", "Total reads", ""),
("", f(m["reads"]), ""),
("R1 too short", "", "R2 too short"),
(d(m["R1_tooshort"],m["reads"]), "", d(m["R2_tooshort"],m["reads"])),
("R1 degen UMI", "", "R2 degen UMI"),
(d(m["R1_homopolymer_UMI"],m["reads"]), "", d(m["R2_homopolymer_UMI"],m["reads"])),
("R1 LQ UMI" , "", "R2 LQ UMI"),
(d(m["R1_N_UMI"],m["reads"]), "", d(m["R2_N_UMI"],m["reads"])),
("R1 no UP", "", "R2 no UP"),
(d(m["R1_no_UP"],m["reads"]), "", d(m["R2_no_UP"],m["reads"])),
("", "Filtered reads", ""),
("", d(m["reads_filtered"], m["reads"]), ""),
("", "Sequencing saturation", ""),
("", s(round((1 - (m["umis_filtered"] / m["reads_filtered"]))*100, digits=1))*"%", ""),
("", "Filtered UMIs", ""),
("", f(m["umis_filtered"]), ""),
("R1 UMI cutoff: ", "", "R2 UMI cutoff"),
(f(m["R1_umicutoff"]), "", f(m["R2_umicutoff"])),
("R1 Barcodes", "R2:R1 ratio", "R2 Barcodes"),
(f(m["R1_barcodes"]), s(round(m["R2_barcodes"]/m["R1_barcodes"],digits=2)), f(m["R2_barcodes"])),
("R1 exact matches", "", "R2 exact matches"),
(d(mm["R1_exact"],m["umis_filtered"]), "", d(mm["R2_exact"],m["umis_filtered"])),
("R1 none matches", "", "R2 none matches"),
(d(mm["R1_none"],m["umis_filtered"]), "", d(mm["R2_none"],m["umis_filtered"])),
("", "Matched UMIs", ""),
("", d(m["umis_matched"], m["umis_filtered"]), ""),
]
p = plot(xlim=(0, 4), ylim=(0, 26+1), framestyle=:none, size=(7*100, 8*100),
         legend=false, xticks=:none, yticks=:none)
for (i, (str1, str2, str3)) in enumerate(data)
    annotate!(p, 1, 26 - i + 1, text(str1, :center, 12))
    annotate!(p, 2, 26 - i + 1, text(str2, :center, 12))
    annotate!(p, 3, 26 - i + 1, text(str3, :center, 12))
end
hline!(p, [12.5], linestyle = :solid, color = :black)
savefig(p, joinpath(out_path, "metadata.pdf"))

merge_pdfs([joinpath(out_path,"elbows.pdf"), joinpath(out_path,"metadata.pdf")], joinpath(out_path,"QC.pdf"), cleanup=true)

@assert isempty(intersect(keys(metadata), keys(matching_metadata)))
meta_df = DataFrame([Dict(:key => k, :value => v) for (k,v) in merge(metadata, matching_metadata)])
sort!(meta_df, :key) ; meta_df = select(meta_df, :key, :value)
CSV.write(joinpath(out_path,"metadata.csv"), meta_df, writeheader=false)

# Write matrix and barcodes
m1 = sort(collect(Set(df.sb1_i)))
m2 = sort(collect(Set(df.sb2_i)))
dict1 = Dict{UInt32, UInt32}(value => index for (index, value) in enumerate(m1))
dict2 = Dict{UInt32, UInt32}(value => index for (index, value) in enumerate(m2))
sb1_whitelist_short = sb1_whitelist[m1]
sb2_whitelist_short = sb2_whitelist[m2]
sb1_new = [dict1[k] for k in df.sb1_i]
sb2_new = [dict2[k] for k in df.sb2_i]
@assert sb1_whitelist[df.sb1_i] == sb1_whitelist_short[sb1_new]
@assert sb2_whitelist[df.sb2_i] == sb2_whitelist_short[sb2_new]
df.sb1_i = sb1_new
df.sb2_i = sb2_new

open(GzipCompressorStream, joinpath(out_path,"sb1.txt.gz"), "w") do file
    for line in sb1_whitelist_short
        write(file, line * "\n")
    end
end
open(GzipCompressorStream, joinpath(out_path,"sb2.txt.gz"), "w") do file
    for line in sb2_whitelist_short
        write(file, line * "\n")
    end
end
open(GzipCompressorStream, joinpath(out_path,"matrix.csv.gz"), "w") do file
    CSV.write(file, df, writeheader=false)
end

println("done!") ; flush(stdout) ; GC.gc()

@assert all(f -> isfile(joinpath(out_path, f)), ["matrix.csv.gz", "sb1.txt.gz", "sb2.txt.gz", "metadata.csv", "QC.pdf"])
