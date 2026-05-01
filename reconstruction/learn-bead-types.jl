using FASTX
using CodecZlib
using DataFrames

# Determine the R1 bead type
function learn_R1type(R1)
    iter = R1 |> open |> GzipDecompressorStream |> FASTQ.Reader
    counts = Dict("V10"=>0, "V17"=>0, "V19"=>0)
    for (i, record) in enumerate(iter)
        i > 100000 ? break : nothing
        seq = FASTQ.sequence(record)
        length(seq) < 36 ? continue : nothing
        counts["V19"] += get_V19(seq)[3] == UP1[1:10]
        length(seq) < 42 ? continue : nothing
        counts["V10"] += get_V10(seq)[3] == UP1
        length(seq) < 44 ? continue : nothing
        counts["V17"] += get_V17(seq)[3] == UP1
    end
    counts["V19"] -= counts["V17"]
    myid() == 1 && println(counts)
    (findmax(counts)[1] < 100000 * 0.01) && error("Unrecognized R1 bead structure for $R1")
    return(findmax(counts)[2])
end

# Determine the R2 bead type
function learn_R2type(R2)
    iter = R2 |> open |> GzipDecompressorStream |> FASTQ.Reader
    counts = Dict("V15"=>0, "V16"=>0)
    for (i, record) in enumerate(iter)
        i > 100000 ? break : nothing
        seq = FASTQ.sequence(record)
        length(seq) < 34 ? continue : nothing
        counts["V15"] += get_V15(seq)[3] == UP2
        length(seq) < 36 ? continue : nothing
        counts["V16"] += get_V16(seq)[3] == UP2
    end
    myid() == 1 && println(counts)
    (findmax(counts)[1] < 100000 * 0.01) && error("Unrecognized R2 bead structure for $R2")
    return(findmax(counts)[2])
end

function learn_bead_types(
    R1s::Vector{String},
    R2s::Vector{String},
)::Tuple{String, String}

    ################################################################################

    R1_types = [learn_R1type(R1) for R1 in R1s]
    if all(x -> x == "V10", R1_types)
        bead1_type = "V10" # const
    elseif all(x -> x == "V17", R1_types)
        bead1_type = "V17" # const
    elseif all(x -> x == "V19", R1_types)
        bead1_type = "V19" # const
    else
        error("The R1 bead type is not consistent ($R1_types)")
    end
    myid() == 1 && println("R1 bead type: $bead1_type")

    R2_types = [learn_R2type(R2) for R2 in R2s]
    if all(x -> x == "V15", R2_types)
        bead2_type = "V15" # const
    elseif all(x -> x == "V16", R2_types)
        bead2_type = "V16" # const
    else
        error("The R2 bead type is not consistent ($R2_types)")
    end
    myid() == 1 && println("R2 bead type: $bead2_type")

    return bead1_type, bead2_type
end
