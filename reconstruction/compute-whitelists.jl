using Plots
using StatsBase
using DataFrames
using Distributions: pdf, Exponential

function compute_cutoffs(
    tab1::Dict{UInt64,Int64},
    tab2::Dict{UInt64,Int64},
    metadata::Dict{String,Int64},
    R1_barcodes::Int64,
    R2_barcodes::Int64,
    out_path::String,
)::Cutoffs

    # Helper method
    function remove_intermediate(x, y)
        m = (y .!= vcat(y[2:end], NaN)) .| (y .!= vcat(NaN, y[1:end-1]))
        x = x[m] ; y = y[m]
        return(x, y)
    end

    # Automatic cutoff
    # Use the elbow plot to determine which beads to use as our whitelist
    #   The cutoff is auto-detected using the steepest part of the curve
    #   To make finding it more consistent, set a reasonable min/max UMI cutoff
    #   uc (umi cutoff) is the steepest part of the curve between min_uc and max_uc
    function determine_umi_cutoff(y)
        sort!(y, rev=true)
        x = 1:length(y)
        x, y = remove_intermediate(x, y)

        # find the steepest slope
        lx = log10.(x) ; ly = log10.(y)
        dydx = (ly[1:end-2] - ly[3:end]) ./ (lx[1:end-2] - lx[3:end])
        min_uc = 10 ; max_uc = 1000 ; m = log10(min_uc) .<= ly[2:end-1] .<= log10(max_uc)
        min_index = findall(m)[argmin(dydx[m])] + 1 + 2

        uc = round(Int64, 10^ly[min_index])
        return uc
    end

    uc1_auto = determine_umi_cutoff(tab1 |> values |> collect) # const
    uc2_auto = determine_umi_cutoff(tab2 |> values |> collect) # const
    bc1_auto = count(e -> e >= uc1_auto, tab1 |> values |> collect) # const
    bc2_auto = count(e -> e >= uc2_auto, tab2 |> values |> collect) # const

    # Manual cutoff
    # R1_barcodes: provided at command-line
    # R2_barcodes: provided at command-line
    function bc_to_uc(bc, table)
        kv = DataFrame(keys = collect(keys(table)), values = collect(values(table)))
        sort!(kv, :keys, rev=true)
        kv[!, :cumsum] = cumsum(kv[!, :values])
        uc = kv[!, :keys][findfirst(x -> x >= bc, kv[!, :cumsum])]
        return uc
    end

    uc1_manual = R1_barcodes > 0 ? bc_to_uc(R1_barcodes, tab1 |> values |> countmap) : R1_barcodes # const
    uc2_manual = R2_barcodes > 0 ? bc_to_uc(R2_barcodes, tab2 |> values |> countmap) : R2_barcodes # const
    bc1_manual = R1_barcodes > 0 ? count(e -> e >= uc1_manual, tab1 |> values |> collect) : R1_barcodes # const
    bc2_manual = R2_barcodes > 0 ? count(e -> e >= uc2_manual, tab2 |> values |> collect) : R2_barcodes # const

    # Use manual cutoff if > 0
    if R1_barcodes > 0
        bc1 = bc1_manual # const
        uc1 = uc1_manual # const
    else
        bc1 = bc1_auto # const
        uc1 = uc1_auto # const
    end
    if R2_barcodes > 0
        bc2 = bc2_manual # const
        uc2 = uc2_manual # const
    else
        bc2 = bc2_auto # const
        uc2 = uc2_auto # const
    end

    # Plots
    function umi_density_plot(table, uc_auto, uc_manual, R)
        x = collect(keys(table))
        y = collect(values(table))
        perm = sortperm(x)
        x = x[perm]
        y = y[perm]

        # Compute the KDE
        lx_s = 0:0.001:ceil(maximum(log10.(x)), digits=3)
        ly_s = []
        for lx_ in lx_s
            weights = [pdf(Exponential(0.05), abs(lx_ - lx)) for lx in log10.(x)]
            kde = sum(log10.(y) .* weights) / sum(weights)
            push!(ly_s, kde)
        end

        # Create a density plot
        p = plot(x, y, seriestype = :scatter, xscale = :log10, yscale = :log10,
                 xlabel = "Number of UMI", ylabel = "Frequency",
                 markersize = 3, markerstrokewidth = 0.1,
                 title = "$R UMI Count Distribution", label = "Barcodes",
                 titlefont=10, guidefont=8, legendfontsize=6, xlabelfontsize=8, ylabelfontsize=8)
        plot!(p, (10).^lx_s, (10).^ly_s, seriestype = :line, label="KDE")
        vline!(p, [uc_auto], linestyle = :dash, color = :red, label = "UMI cutoff (auto)")
        uc_manual > 0 && vline!(p, [uc_manual], linestyle = :dash, color = :purple, label = "UMI cutoff (manual)")
        xticks!(p, [10^i for i in 0:ceil(log10(maximum(x)))])
        yticks!(p, [10^i for i in 0:ceil(log10(maximum(y)))])
        return p
    end

    p1 = umi_density_plot(tab1 |> values |> countmap, uc1_auto, uc1_manual, "R1")
    p3 = umi_density_plot(tab2 |> values |> countmap, uc2_auto, uc2_manual, "R2")

    function elbow_plot(y, uc_auto, uc_manual, bc_auto, bc_manual, R)
        sort!(y, rev=true)
        x = 1:length(y)

        xp, yp = remove_intermediate(x, y)
        p = plot(xp, yp, seriestype = :line, xscale = :log10, yscale = :log10,
             xlabel = "$R Spatial Barcode Rank", ylabel = "UMI Count",
             title = "$R Spatial Barcode Elbow Plot", label = "Barcodes",
             titlefont=10, guidefont=8, legendfontsize=6, xlabelfontsize=8, ylabelfontsize=8)
        hline!(p, [uc_auto], linestyle = :dash, color = :red, label = "UMI cutoff (auto)")
        uc_manual > 0 && hline!(p, [uc_manual], linestyle = :dash, color = :purple, label = "UMI cutoff (manual)")
        vline!(p, [bc_auto], linestyle = :dash, color = :green, label = "SB cutoff (auto)")
        bc_manual > 0 && vline!(p, [bc_manual], linestyle = :dash, color = :brown, label = "SB cutoff (manual)")
        xticks!(p, [10^i for i in 0:ceil(log10(maximum(xp)))])
        yticks!(p, [10^i for i in 0:ceil(log10(maximum(yp)))])
        return p
    end

    p2 = elbow_plot(tab1 |> values |> collect, uc1_auto, uc1_manual, bc1_auto, bc1_manual, "R1")
    p4 = elbow_plot(tab2 |> values |> collect, uc2_auto, uc2_manual, bc2_auto, bc2_manual, "R2")

    p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
    savefig(p, joinpath(out_path, "elbows.pdf"))

    metadata["R1_umicutoff_auto"] = uc1_auto
    metadata["R2_umicutoff_auto"] = uc2_auto
    metadata["R1_barcodes_auto"] = bc1_auto
    metadata["R2_barcodes_auto"] = bc2_auto

    metadata["R1_umicutoff_manual"] = uc1_manual
    metadata["R2_umicutoff_manual"] = uc2_manual
    metadata["R1_barcodes_manual"] = bc1_manual
    metadata["R2_barcodes_manual"] = bc2_manual

    metadata["R1_umicutoff"] = uc1
    metadata["R2_umicutoff"] = uc2
    metadata["R1_barcodes"] = bc1
    metadata["R2_barcodes"] = bc2

    return Cutoffs(uc1, bc1, uc2, bc2)
end

function compute_whitelists(
    df::DataFrame,
    metadata::Dict{String,Int64},
    R1_barcodes::Int64,
    R2_barcodes::Int64,
    out_path::String,
)::Tuple{Dict{UInt64,Int64}, Dict{UInt64,Int64}, Cutoffs}
    print_start("Computing barcode whitelist... ")

    # Count the number of times each barcode appears
    tab1 = countmap(df[!,:sb1_i]) # const
    tab2 = countmap(df[!,:sb2_i]) # const

    cutoffs = compute_cutoffs(tab1, tab2, metadata, R1_barcodes, R2_barcodes, out_path)

    println_done()
    return tab1, tab2, cutoffs
end

function compute_whitelists(
    tab1::Dict{UInt64,Int64},
    tab2::Dict{UInt64,Int64},
    metadata::Dict{String,Int64},
    R1_barcodes::Int64,
    R2_barcodes::Int64,
    out_path::String,
)::Cutoffs
    print_start("Computing barcode whitelist... ")

    cutoffs = compute_cutoffs(tab1, tab2, metadata, R1_barcodes, R2_barcodes, out_path)

    println_done()
    return cutoffs
end
