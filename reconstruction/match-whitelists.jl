using DataFrames

function make_whitelist(
    tab1::Dict{UInt64,Int64},
    tab2::Dict{UInt64,Int64},
    cutoffs::Cutoffs,
)::Tuple{Set{UInt64}, Set{UInt64}}
    uc1 = cutoffs.uc1
    bc1 = cutoffs.bc1
    uc2 = cutoffs.uc2
    bc2 = cutoffs.bc2
    wl1 = Set{UInt64}([k for (k, v) in tab1 if v >= uc1]) ; @assert length(wl1) == bc1
    wl2 = Set{UInt64}([k for (k, v) in tab2 if v >= uc2]) ; @assert length(wl2) == bc2
    return wl1, wl2
end

function match_barcode(df, metadata, wl1, wl2)
    m1 = [s1 in wl1 for s1 in df.sb1_i]
    m2 = [s2 in wl2 for s2 in df.sb2_i]

    metadata["R1_exact"] = sum(m1)
    metadata["R2_exact"] = sum(m2)

    df.keep = m1 .& m2
    filter!(:keep => identity, df)
    select!(df, Not(:keep))

    return nothing
end

function match_whitelists(
    df::DataFrame,
    metadata::Dict{String,Int64},
    tab1::Dict{UInt64,Int64},
    tab2::Dict{UInt64,Int64},
    cutoffs::Cutoffs,
)::Nothing

    print_start("Matching to barcode whitelist... ")

    wl1, wl2 = make_whitelist(tab1, tab2, cutoffs)
    match_barcode(df, metadata, wl1, wl2) # this function modifies in-place
    metadata["umis_exact"] = nrow(df)

    println_done()
    return nothing
end

function save_whitelists(
    tab1::Dict{UInt64,Int64},
    tab2::Dict{UInt64,Int64},
    cutoffs::Cutoffs,
    out_path::String,
)::Nothing

    print_start("Matching to barcode whitelist... ")

    wl1, wl2 = make_whitelist(tab1, tab2, cutoffs)

    println_done()

    print_start("Writing barcode whitelist... ")

    write_barcode_set(wl1, joinpath(out_path, "wl1.txt.gz"))
    write_barcode_set(wl2, joinpath(out_path, "wl2.txt.gz"))

    println_done()

    return nothing
end

function match_whitelists(
    df::DataFrame,
    metadata::Dict{String,Int64},
    in_path::String,
)::Nothing

    print_start("Reading barcode whitelist... ")

    wl1 = read_barcode_set(joinpath(in_path, "wl1.txt.gz"))
    wl2 = read_barcode_set(joinpath(in_path, "wl2.txt.gz"))

    println_done()

    print_start("Matching to barcode whitelist... ")

    match_barcode(df, metadata, wl1, wl2) # this function modifies in-place
    metadata["umis_exact"] = nrow(df)

    println_done()
    return nothing
end
