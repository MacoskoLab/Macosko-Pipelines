using ArgParse
using Distributed

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

################################################################################

include(joinpath(@__DIR__, "logging.jl"))

include(joinpath(@__DIR__, "find-fastqs.jl"))

const R1s, R2s = find_fastqs(fastq_path, regex, out_path)

################################################################################

# Create a worker for each FASTQ pair
addprocs(length(R1s))

@everywhere include(joinpath(@__DIR__, "bead-info.jl"))

include(joinpath(@__DIR__, "learn-bead-types.jl"))

const bead1_type, bead2_type = learn_bead_types(R1s, R2s)

include(joinpath(@__DIR__, "read-fastqs.jl"))

const metadatas = process_results(out_path, prob, bead1_type, bead2_type, R1s, R2s)

rmprocs(workers())

const df, metadata = combine_results(out_path, metadatas)

################################################################################

const bead1_info = bead1_type_to_info(bead1_type)
const bead2_info = bead2_type_to_info(bead2_type)

include(joinpath(@__DIR__, "data-frame.jl"))

include(joinpath(@__DIR__, "count-reads.jl"))

count_reads(df, metadata, bead1_info, bead2_info, out_path)

################################################################################

# # Save results of fastq parsing
# h5open(joinpath(out_path, "reads.h5"), "w") do file
#     file["sb1_2bit", compress=1] = df[!, :sb1_i]
#     file["umi1_2bit", compress=1] = df[!, :umi1_i]
#     file["sb2_2bit", compress=1] = df[!, :sb2_i]
#     file["umi2_2bit", compress=1] = df[!, :umi2_i]
#     file["reads", compress=1] = df[!, :reads]
# end

# # Load previous fastq parsing results
# df = h5open(joinpath(out_path, "reads.h5"), "r") do file 
#     DataFrame(sb1_i = read(file["sb1_2bit"]),
#               umi1_i = read(file["umi1_2bit"]),
#               sb2_i = read(file["sb2_2bit"]),
#               umi2_i = read(file["umi2_2bit"]),
#               reads = read(file["reads"]))
# end
# metadata = Dict(String(row[1]) => parse(Int, row[2]) for row in CSV.Rows(joinpath(out_path, "metadata.csv"), header=false))

################################################################################

include(joinpath(@__DIR__, "cutoffs.jl"))

include(joinpath(@__DIR__, "compute-whitelists.jl"))

const tab1, tab2, cutoffs = compute_whitelists(df, metadata, R1_barcodes, R2_barcodes, out_path)

################################################################################

include(joinpath(@__DIR__, "barcode-set.jl"))

include(joinpath(@__DIR__, "match-whitelists.jl"))

match_whitelists(df, metadata, tab1, tab2, cutoffs)

################################################################################

include(joinpath(@__DIR__, "remove-chimeras.jl"))

compute_chimeric1(df)
compute_chimeric2(df)
remove_chimeras(df, metadata)

################################################################################

include(joinpath(@__DIR__, "count-umis.jl"))

count_umis(df, out_path)

################################################################################

include(joinpath(@__DIR__, "filter-connections.jl"))

filter_connections(df, metadata, out_path)

################################################################################

include(joinpath(@__DIR__, "write-outputs.jl"))

write_outputs(df, metadata, R1_barcodes, R2_barcodes, bead1_info, bead2_info, prob, out_path)
