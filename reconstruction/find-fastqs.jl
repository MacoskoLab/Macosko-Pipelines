using Plots

function find_fastqs(
    fastq_path::String,
    regex::Regex,
    out_path::String,
)::Tuple{Vector{String}, Vector{String}}
    # Load the FASTQ paths
    fastqs = readdir(fastq_path, join=true)
    fastqs = filter(fastq -> endswith(fastq, ".fastq.gz"), fastqs)
    fastqs = filter(fastq -> occursin(regex, fastq), fastqs)
    @assert length(fastqs) >= 2 "ERROR: No FASTQ pairs found"
    R1s = filter(s -> occursin("_R1_001", s), fastqs) ; println("R1s: ", basename.(R1s)) # const
    R2s = filter(s -> occursin("_R2_001", s), fastqs) ; println("R2s: ", basename.(R2s)) # const
    @assert length(R1s) > 0 && length(R2s) > 0 "ERROR: No FASTQ pairs found"
    @assert length(R1s) == length(R2s) "ERROR: R1s and R2s are not all paired"
    @assert [replace(R1, "_R1_001"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_001"=>"", count=1) for R2 in R2s]
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
    return R1s, R2s
end
