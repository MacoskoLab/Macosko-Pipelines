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

const in_path = args["in_path"]
println("Input path: "*in_path)
@assert isdir(in_path) "Input path not found"
@assert !isempty(readdir(in_path)) "Input path is empty"

const out_path = args["out_path"]
println("Output path: "*out_path)
Base.Filesystem.mkpath(out_path)
@assert isdir(out_path) "Output path could not be created"

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

include(joinpath(@__DIR__, "metadata.jl"))

const metadata = read_metadata(in_path)

################################################################################

include(joinpath(@__DIR__, "bead-info.jl"))

include(joinpath(@__DIR__, "data-frame.jl"))

include(joinpath(@__DIR__, "count-reads.jl"))

const tab1, tab2 = read_rpsbs(in_path)

################################################################################

include(joinpath(@__DIR__, "cutoffs.jl"))

include(joinpath(@__DIR__, "compute-whitelists.jl"))

const cutoffs = compute_whitelists(tab1, tab2, metadata, R1_barcodes, R2_barcodes, out_path)

################################################################################

write_metadata(out_path, metadata)

################################################################################

include(joinpath(@__DIR__, "barcode-set.jl"))

include(joinpath(@__DIR__, "match-whitelists.jl"))

save_whitelists(tab1, tab2, cutoffs, out_path)
