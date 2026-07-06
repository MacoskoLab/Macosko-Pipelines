using ArgParse
using DataFrames
using Distributed

# Count the reads only for the specified prefixes.
# Writes a plot showing the file paths and the prefix.

# Read the command-line arguments
function get_args()
    s = ArgParseSettings()

    # Positional arguments
    @add_arg_table s begin
        "in_path"
        help = "Input directory"
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

        "--R1_filter", "-f"
        help = "Include only R1 reads that start with this prefix"
        arg_type = String
        default = ""

        "--R2_filter", "-g"
        help = "Include only R2 reads that start with this prefix"
        arg_type = String
        default = ""
    end

    return parse_args(ARGS, s)
end

# Load the command-line arguments
args = get_args()

const in_path = args["in_path"]
println("Input path: "*in_path)
@assert isdir(in_path) "Input path not found"
@assert !isempty(readdir(in_path)) "Input path is empty"

const out_path = args["out_path"]
println("Output path: "*out_path)
Base.Filesystem.mkpath(out_path)
@assert isdir(out_path) "Output path could not be created"

const prob = args["downsampling_level"]
@assert 0 < prob <= 1 "Invalid downsampling level $prob"
if prob < 1
    println("Downsampling level: $prob")
end

const R1_filter = args["R1_filter"]
const R2_filter = args["R2_filter"]

if length(R1_filter) > 0 && length(R2_filter) > 0
    const file_suffix = ".R1_$R1_filter.R2_$R2_filter"
elseif length(R1_filter) > 0
    const file_suffix = ".R1_$R1_filter"
elseif length(R2_filter) > 0
    const file_suffix = ".R2_$R2_filter"
else
    const file_suffix = ""
end

if length(R1_filter) > 0 || length(R2_filter) > 0
    println("R1 filter: $R1_filter*")
    println("R2 filter: $R2_filter*")
    println("File suffix: $file_suffix")
end

println("Threads: $(Threads.nthreads())\n")

################################################################################

include(joinpath(@__DIR__, "logging.jl"))

include(joinpath(@__DIR__, "find-fastqs.jl"))

const R1s, R2s = load_fastqs(in_path)

################################################################################

# Create a worker for each FASTQ pair
addprocs(length(R1s))

# After the workers are created, we can include the necessary files on all workers.
# Let each worker and the main code know about bead info.
@everywhere include(joinpath(@__DIR__, "bead-info.jl"))

include(joinpath(@__DIR__, "metadata.jl"))

include(joinpath(@__DIR__, "learn-bead-types.jl"))

const bead1_type, bead2_type = metadata_to_bead_types(read_metadata(in_path))

include(joinpath(@__DIR__, "data-frame.jl"))

include(joinpath(@__DIR__, "read-fastqs.jl"))

# Read the fastqs on the workers
const metadatas = process_results(out_path, prob, bead1_type, bead2_type, R1s, R2s, R1_filter, R2_filter)

# We're done with the distributed workers so remove them.
rmprocs(workers())

const df, metadata = combine_results(out_path, metadatas)

################################################################################

const bead1_info = bead1_type_to_info(bead1_type)
const bead2_info = bead2_type_to_info(bead2_type)

include(joinpath(@__DIR__, "count-reads.jl"))

count_reads(df, metadata, bead1_info, bead2_info, out_path, file_suffix)

################################################################################

write_rpsbs(df, out_path, file_suffix)

################################################################################

include(joinpath(@__DIR__, "metadata.jl"))

write_metadata(out_path, metadata, file_suffix)

################################################################################

write_reads_df(df, out_path, file_suffix)
