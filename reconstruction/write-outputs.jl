using CSV
using Plots
using PDFmerger
using CodecZlib
using DataFrames
using DelimitedFiles

function write_outputs(
    df::DataFrame,
    metadata::Dict{String,Int64},
    R1_barcodes::Int64,
    R2_barcodes::Int64,
    bead1_info::BeadInfo,
    bead2_info::BeadInfo,
    prob::Float64,
    out_path::String,
)::Nothing
    bead1_type = bead1_info.bead_type
    bead2_type = bead2_info.bead_type
    decode_sb1 = bead1_info.decode_sb
    decode_sb2 = bead2_info.decode_sb

    print("Writing output... ") ; flush(stdout)

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

    merge_pdfs([joinpath(out_path,"elbows.pdf"),
                joinpath(out_path,"metadata.pdf"),
                joinpath(out_path,"SNR.pdf"),
                joinpath(out_path,"connection_filter.pdf"),
                joinpath(out_path,"histograms.pdf"),
                joinpath(out_path,"filepaths.pdf")],
                joinpath(out_path,"QC.pdf"), cleanup=true)

    # Factorize the barcode indexes
    uniques1 = sort(collect(Set(df.sb1_i))) # called sb1_i (encoded)
    uniques2 = sort(collect(Set(df.sb2_i))) # called sb2_i (encoded)
    sb1_whitelist = [decode_sb1(sb1_i) for sb1_i in uniques1] # called sb1 (barcodes)
    sb2_whitelist = [decode_sb2(sb2_i) for sb2_i in uniques2] # called sb2 (barcodes)
    dict1 = Dict{UInt64, UInt64}(value => index for (index, value) in enumerate(uniques1))
    dict2 = Dict{UInt64, UInt64}(value => index for (index, value) in enumerate(uniques2))
    df.sb1_i = [dict1[k] for k in df.sb1_i]
    df.sb2_i = [dict2[k] for k in df.sb2_i]
    @assert sort(collect(Set(df.sb1_i))) == collect(1:length(Set(df.sb1_i)))
    @assert sort(collect(Set(df.sb2_i))) == collect(1:length(Set(df.sb2_i)))

    # Save the matrix
    open(GzipCompressorStream, joinpath(out_path, "sb1.txt.gz"), "w") do file
        writedlm(file, sb1_whitelist, "\n")
    end
    open(GzipCompressorStream, joinpath(out_path, "sb2.txt.gz"), "w") do file
        writedlm(file, sb2_whitelist, "\n")
    end
    rename!(df, Dict(:sb1_i => :sb1_index, :sb2_i => :sb2_index, :umi => :umi))
    CSV.write(joinpath(out_path, "matrix.csv.gz"), df, writeheader=true, compress=true)

    @assert all(f -> isfile(joinpath(out_path, f)), ["matrix.csv.gz", "sb1.txt.gz", "sb2.txt.gz", "QC.pdf", "metadata.csv"])

    println("done") ; flush(stdout) ; GC.gc()
    return nothing
end
