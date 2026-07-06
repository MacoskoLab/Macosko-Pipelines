using ArgParse
using Distributed

# Find the fastq paths and write out a plot showing the file paths.
# The plot generated is added to the final QC.pdf.

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

println("Threads: $(Threads.nthreads())\n")

################################################################################

include(joinpath(@__DIR__, "logging.jl"))

include(joinpath(@__DIR__, "find-fastqs.jl"))

const R1s, R2s = find_fastqs(fastq_path, regex, out_path)

################################################################################

include(joinpath(@__DIR__, "bead-info.jl"))

include(joinpath(@__DIR__, "learn-bead-types.jl"))

const bead1_type, bead2_type = learn_bead_types(R1s, R2s)

################################################################################

const metadata = bead_types_to_metadata(bead1_type, bead2_type)

################################################################################

include(joinpath(@__DIR__, "metadata.jl"))

write_metadata(out_path, metadata)

################################################################################

save_fastqs(R1s, R2s, out_path)
