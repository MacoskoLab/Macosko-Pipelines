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

const metadata = read_metadatas(in_path)

################################################################################

include(joinpath(@__DIR__, "data-frame.jl"))

include(joinpath(@__DIR__, "count-umis.jl"))

const df = plot_umis(in_path, out_path)

################################################################################

include(joinpath(@__DIR__, "filter-connections.jl"))

filter_connections(df, metadata, out_path)

################################################################################

include(joinpath(@__DIR__, "bead-info.jl"))

const bead1_info, bead2_info = metadata_to_bead_infos(metadata)

################################################################################

include(joinpath(@__DIR__, "write-outputs.jl"))

# If the input path is not the same as the output path, copy the PDFs to the output path.
if in_path != out_path
    copy_pdfs(in_path, out_path)
end

################################################################################

write_outputs(df, metadata, R1_barcodes, R2_barcodes, bead1_info, bead2_info, prob, out_path)
