using ArgParse
using Distributed

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
        "--R1_filter", "-f"
        help = "Only compute chimeras for R1 reads that start with this prefix"
        arg_type = String
        default = ""

        "--R2_filter", "-g"
        help = "Only compute chimeras for R2 reads that start with this prefix"
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

const R1_filter = args["R1_filter"]
const R2_filter = args["R2_filter"]

# Both filters are not allowed to be set, otherwise we might miss chimeras due to the splitting of the data frame.
if length(R1_filter) > 0 && length(R2_filter) > 0
    error("ERROR: Both R1 and R2 filters cannot be set at the same time")
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

include(joinpath(@__DIR__, "bead-info.jl"))

include(joinpath(@__DIR__, "data-frame.jl"))

include(joinpath(@__DIR__, "count-reads.jl"))

const df = read_reads_dfs(in_path, file_suffix)

################################################################################

include(joinpath(@__DIR__, "remove-chimeras.jl"))

# If the R2 filter is empty, then we need to compute the chimeras for the R1 reads
if length(R2_filter) == 0
    compute_chimeric1(df)
    write_chimeric1(df, out_path, file_suffix)
end

# If the R1 filter is empty, then we need to compute the chimeras for the R2 reads
if length(R1_filter) == 0
    compute_chimeric2(df)
    write_chimeric2(df, out_path, file_suffix)
end
