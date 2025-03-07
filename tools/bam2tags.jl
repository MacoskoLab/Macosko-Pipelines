using XAM
using ArgParse
using CodecZlib

# Read the command-line arguments
function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "bam_path"
        help = "Input BAM path/name (.bam)"
        arg_type = String
        required = true

        "out_path"
        help = "Output file path/name (.tsv.gz)"
        arg_type = String
        required = true
    end

    return parse_args(ARGS, s)
end

# Load the command-line arguments
args = get_args()

const bampath = args["bam_path"]
@assert isfile(bampath) "ERROR: BAM could not be found"

const outpath = args["out_path"]
Base.Filesystem.mkpath(dirname(outpath))
dirname(outpath) != "" && @assert isdir(dirname(outpath)) "ERROR: Output path could not be created"

# Get all the unique tags
taghash(s) = Int(s[1]-48) + Int(s[2]-48)*75 + 1
tagorder = ["NH","HI","NM","MD","AS","nM","jM","jI",
            "uT","RG",
            "CB","CR","CY","UB","UR","UY","TR",
            "TX","AN","GX","GN","MM","RE","pa","pr","ts","xf",
            "fb","fr","fq","fx"]
println("Getting all unique tags...")
cmd = pipeline(`samtools view "$bampath"`,`cut -f 12-`,`tr '\t' '\n'`,`cut -d':' -f 1`,`awk '!a[$0]++'`)
uniquetags = [string(s) for s in split(strip(read(cmd, String)), '\n')]
tags = [[tag for tag in tagorder if tag in uniquetags]; [tag for tag in uniquetags if tag ∉ tagorder]]
const indexes = zeros(Int, taghash("zz"))
for (i, tag) in enumerate(tags)
    indexes[taghash(tag)] = i
end
const ntags = length(tags)

# Parse the BAM
println("Beginning loop...")
@time begin
    out = GzipCompressorStream(open(outpath, "w"))
    write(out, join(tags, '\t') * '\n')
    reader = open(BAM.Reader, bampath)

    for record in reader
        row = fill("", ntags)
        for tag in BAM.auxdata(record)
            row[indexes[taghash(tag[1])]] = string(tag[2])
        end
        write(out, join(row, '\t') * '\n')
    end
    
    close(reader)
    close(out)
end

println("Done!")
