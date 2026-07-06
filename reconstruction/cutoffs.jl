# UMI cutoffs and barcode counts that made the cutoff.
# The actual list of barcodes is saved in a separate file.

struct Cutoffs
    uc1::Int64
    bc1::Int64
    uc2::Int64
    bc2::Int64
end

function metadata_to_cutoffs(metadata::Dict{String,Int64})::Cutoffs
    uc1 = metadata["R1_umicutoff"]
    bc1 = metadata["R1_barcodes"]
    uc2 = metadata["R2_umicutoff"]
    bc2 = metadata["R2_barcodes"]
    return Cutoffs(uc1, bc1, uc2, bc2)
end
