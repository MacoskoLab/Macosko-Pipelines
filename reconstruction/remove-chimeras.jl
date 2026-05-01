using DataFrames

function remove_chimeras(
    df::DataFrame,
    metadata::Dict{String,Int64},
)::Nothing
    print("Removing chimeras... ") ; flush(stdout)

    #function remove_chimeras(df, metadata)
        if nrow(df) == 0
            metadata["R1_chimeric"] = 0
            metadata["R2_chimeric"] = 0
            return
        end

        sort!(df, [:sb1_i, :umi1_i, :reads], rev = [false, false, true])
        before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb1_i,:umi1_i]]))
        after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb1_i,:umi1_i,:reads]]), false)

        df.chimeric1 = before_same .| after_same

        sort!(df, [:sb2_i, :umi2_i, :reads], rev = [false, false, true])
        before_same = vcat(false, reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb2_i,:umi2_i]]))
        after_same = vcat(reduce(.&, [df[2:end,c] .== df[1:end-1,c] for c in [:sb2_i,:umi2_i,:reads]]), false)
        df.chimeric2 = before_same .| after_same

        metadata["R1_chimeric"] = sum(df.chimeric1)
        metadata["R2_chimeric"] = sum(df.chimeric2)

        subset!(df, :chimeric1 => x -> .!x, :chimeric2 => x -> .!x)
        select!(df, Not([:chimeric1, :chimeric2]))
    #    nothing
    #end
    #remove_chimeras(df, metadata) # this function modifies in-place
    metadata["umis_chimeric"] = metadata["umis_exact"] - nrow(df)

    println("done") ; flush(stdout) ; GC.gc()
    return nothing
end
