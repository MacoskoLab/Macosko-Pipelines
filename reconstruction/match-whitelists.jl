using DataFrames

function match_whitelists(
    df::DataFrame,
    metadata::Dict{String,Int64},
    tab1::Dict{UInt64,Int64},
    tab2::Dict{UInt64,Int64},
    cutoffs::Cutoffs,
)::Nothing

    uc1 = cutoffs.uc1
    bc1 = cutoffs.bc1
    uc2 = cutoffs.uc2
    bc2 = cutoffs.bc2

    print("Matching to barcode whitelist... ") ; flush(stdout)

    function match_barcode(df, metadata)
        wl1 = Set{UInt64}([k for (k, v) in tab1 if v >= uc1]) ; @assert length(wl1) == bc1
        wl2 = Set{UInt64}([k for (k, v) in tab2 if v >= uc2]) ; @assert length(wl2) == bc2

        m1 = [s1 in wl1 for s1 in df.sb1_i]
        m2 = [s2 in wl2 for s2 in df.sb2_i]

        metadata["R1_exact"] = sum(m1)
        metadata["R2_exact"] = sum(m2)

        df.keep = m1 .& m2
        filter!(:keep => identity, df)
        select!(df, Not(:keep))

        nothing
    end
    match_barcode(df, metadata) # this function modifies in-place
    metadata["umis_exact"] = nrow(df)

    println("done") ; flush(stdout) ; GC.gc()
    return nothing
end
