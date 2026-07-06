using StringViews

struct BeadInfo
    bead_type::String
    R_len::Int
    get_R::Function
    encode_sb::Function
    decode_sb::Function
end

# Read structure methods
const SeqView = StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}}
@inline function get_V10(seq::SeqView)
    @inbounds sb_1 = seq[1:8]
    @inbounds up = seq[9:26]
    @inbounds sb_2 = seq[27:33]
    @inbounds umi = seq[34:42]
    return sb_1, sb_2, up, umi
end
@inline function get_V17(seq::SeqView)
    @inbounds sb_1 = seq[1:9]
    @inbounds up = seq[10:27]
    @inbounds sb_2 = seq[28:35]
    @inbounds umi = seq[36:44]
    return sb_1, sb_2, up, umi
end
@inline function get_V19(seq::SeqView)
    @inbounds sb_1 = seq[1:9]
    @inbounds up = seq[10:19]
    @inbounds sb_2 = seq[20:27]
    @inbounds umi = seq[28:36]
    return sb_1, sb_2, up, umi
end
@inline function get_V15(seq::SeqView)
    @inbounds sb_1 = seq[1:8]
    @inbounds sb_2 = seq[9:15]
    @inbounds up = seq[16:25]
    @inbounds umi = seq[26:34]
    return sb_1, sb_2, up, umi
end
@inline function get_V16(seq::SeqView)
    @inbounds sb_1 = seq[1:9]
    @inbounds sb_2 = seq[10:17]
    @inbounds up = seq[18:27]
    @inbounds umi = seq[28:36]
    return sb_1, sb_2, up, umi
end
const UP1 = "TCTTCAGCGTTCCCGAGA"
const UP2 = "CTGTTTCCTG"

# String bit-encoding methods
const bases = ['A','C','T','G'] # MUST NOT change this order
const px7 = [convert(UInt32, 4^i) for i in 0:6]
const px8 = [convert(UInt32, 4^i) for i in 0:7]
const px9 = [convert(UInt32, 4^i) for i in 0:8]

@inline function encode_str(str::String)::UInt64 # careful, encodes N as G
    return dot([4^i for i in 0:(length(str)-1)], (codeunits(str) .>> 1) .& 3)
end

@inline function encode_umi(umi::SeqView)::UInt32
    @fastmath @inbounds b = dot(px9, (codeunits(umi) .>> 1) .& 3)
    return b
end
@inline function decode_umi(code::UInt32)::String
    @fastmath @inbounds u = [bases[(code >> n) & 3 + 1] for n in 0:2:16]
    return String(u)
end

@inline function encode_15(sb_1::SeqView, sb_2::SeqView)::UInt64
    @fastmath @inbounds b1 = dot(px8, (codeunits(sb_1) .>> 1) .& 3)
    @fastmath @inbounds b2 = dot(px7, (codeunits(sb_2) .>> 1) .& 3)
    return b1 + b2 * 4^8
end
@inline function decode_15(code::UInt64)::String
    @fastmath @inbounds u = [bases[(code >> n) & 3 + 1] for n in 0:2:28]
    return String(u)
end

@inline function encode_17(sb_1::SeqView, sb_2::SeqView)::UInt64
    @fastmath @inbounds b1 = dot(px9, (codeunits(sb_1) .>> 1) .& 3)
    @fastmath @inbounds b2 = dot(px8, (codeunits(sb_2) .>> 1) .& 3)
    return b1 + b2 * 4^9
end
@inline function decode_17(code::UInt64)::String
    @fastmath @inbounds u = [bases[(code >> n) & 3 + 1] for n in 0:2:32]
    return String(u)
end

function bead1_type_to_info(bead1_type::String)::BeadInfo
    if bead1_type == "V10"
        R1_len = 42
        get_R1 = get_V10
        encode_sb1 = encode_15
        decode_sb1 = decode_15
    elseif bead1_type == "V17"
        R1_len = 44
        get_R1 = get_V17
        encode_sb1 = encode_17
        decode_sb1 = decode_17
    elseif bead1_type == "V19"
        R1_len = 36
        get_R1 = get_V19
        encode_sb1 = encode_17
        decode_sb1 = decode_17
    else
        error("Unrecognized R1 bead type: $bead1_type")
    end

    return BeadInfo(bead1_type, R1_len, get_R1, encode_sb1, decode_sb1)
end

function bead2_type_to_info(bead2_type::String)::BeadInfo
    if bead2_type == "V15"
        R2_len = 34
        get_R2 = get_V15
        encode_sb2 = encode_15
        decode_sb2 = decode_15
    elseif bead2_type == "V16"
        R2_len = 36
        get_R2 = get_V16
        encode_sb2 = encode_17
        decode_sb2 = decode_17
    else
        error("Unrecognized R2 bead type: $bead2_type")
    end

    return BeadInfo(bead2_type, R2_len, get_R2, encode_sb2, decode_sb2)
end

function bead_types_to_metadata(bead1_type::String, bead2_type::String)::Dict{String,Int64}
    metadata = Dict{String,Int64}()
    metadata["R1_beadtype"] = parse(Int, join(filter(isdigit, bead1_type)))
    metadata["R2_beadtype"] = parse(Int, join(filter(isdigit, bead2_type)))
    return metadata
end

function metadata_to_bead_types(metadata::Dict{String,Int64})::Tuple{String, String}
    bead1_type = "V"*string(metadata["R1_beadtype"])
    bead2_type = "V"*string(metadata["R2_beadtype"])
    return bead1_type, bead2_type
end

function metadata_to_bead_infos(metadata::Dict{String,Int64})::Tuple{BeadInfo, BeadInfo}
    bead1_type, bead2_type = metadata_to_bead_types(metadata)
    bead1_info = bead1_type_to_info(bead1_type)
    bead2_info = bead2_type_to_info(bead2_type)
    return bead1_info, bead2_info
end
