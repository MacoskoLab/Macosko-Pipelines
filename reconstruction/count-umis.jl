using Plots
using StatsBase
using StatsPlots
using DataFrames

function count_umis(
    df::DataFrame,
    out_path::String,
)::Nothing
    print("Counting UMIs... ") ; flush(stdout)

    function hist_10cap(vec)
        tab = countmap(vec)
        plotdf = DataFrame(value = collect(keys(tab)), count = collect(values(tab)))
        sum10 = sum(filter(:value => v -> v >= 10, plotdf).count)
        filter!(:value => v -> v < 10, plotdf)
        push!(plotdf, (10, sum10))
        p = bar(plotdf.value, plotdf.count, legend = false,
                xticks = (1:10, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"]),
                titlefont = 10, guidefont = 8)
        return p
    end
    p1 = hist_10cap(df.reads)
    xlabel!(p1, "Reads per UMI")
    ylabel!(p1, "Number of filtered UMIs")
    title!(p1, "Read depth")

    function count_umis(df)
        select!(df, [:sb1_i, :sb2_i])
        if nrow(df) == 0
            df.umi = UInt64[]
            return
        end
        sort!(df, [:sb1_i, :sb2_i])
        start = vcat(true, (df.sb1_i[2:end] .!= df.sb1_i[1:end-1]) .| (df.sb2_i[2:end] .!= df.sb2_i[1:end-1]))
        umis = vcat(diff(findall(start)), nrow(df)-findlast(start)+1)

        df.start = start ; filter!(:start => identity, df) ; select!(df, Not(:start))
        df.umi = umis
        nothing
    end
    count_umis(df) # this function modifies in-place

    p2 = hist_10cap(df.umi)
    xlabel!(p2, "UMIs per connection")
    ylabel!(p2, "Number of connections")
    title!(p2, "Connection distribution")

    function sum_topn(v, n)
        return sum(sort(v, rev=true)[1:min(n, length(v))])
    end
    function plot_umi_distributions(df, col::Symbol)
        R = "R"*string(col)[3] # R1 or R2

        if nrow(df) == 0
            p1 = plot(title="$R SNR", titlefont=10, guidefont=8, xticks=(1:4, ["5", "20", "50", "100"]), yticks=[0.0, 0.25, 0.5, 0.75, 1.0], xlabel="Number of top beads", ylabel="%UMIs in top beads")
            p2 = plot(title="$R UMI Distribution", titlefont=10, guidefont=8, xlabel="log10 UMIs", ylabel="log10 connections")
            return p1, p2
        end

        gdf = combine(groupby(df, col), :umi => sum => :umi,
                                        :umi => length => :connections,
                                        :umi => (v->sum_topn(v,5)) => :top5,
                                        :umi => (v->sum_topn(v,20)) => :top20,
                                        :umi => (v->sum_topn(v,50)) => :top50,
                                        :umi => (v->sum_topn(v,100)) => :top100)

        plotdf = vcat(DataFrame(x = 1, y = gdf.top5 ./ gdf.umi),
                      DataFrame(x = 2, y = gdf.top20 ./ gdf.umi),
                      DataFrame(x = 3, y = gdf.top50 ./ gdf.umi),
                      DataFrame(x = 4, y = gdf.top100 ./ gdf.umi))

        # SNR violins
        p1 = @df plotdf begin
            violin(:x, :y, line = 0, fill = (0.3, :blue), legend = false, titlefont = 10, guidefont = 8,
                xticks = ([1, 2, 3, 4], ["5", "20", "50", "100"]), yticks = [0.0, 0.25, 0.5, 0.75, 1.0],
            xlabel = "Number of top beads", ylabel = "%UMIs in top beads", title = "$R SNR")
            boxplot!(:x, :y, line = (1, :black), fill = (0.3, :grey), outliers = false, legend = false)
        end

        # log-umi vs log-connection density
        m = max(log10(maximum(gdf.umi)),log10(maximum(gdf.connections)))
        num_bins = length(unique(gdf.umi)) < 2 ? 2 : :auto # Avoid ERROR: LoadError: InexactError: Int64(NaN)
        xmean = round(log10(mean(gdf.umi)), digits=2)
        xmed = round(log10(median(gdf.umi)), digits=2)
        ymean = round(log10(mean(gdf.connections)), digits=2)
        ymed = round(log10(median(gdf.connections)),digits=2)
        p2 = histogram2d(log10.(gdf.umi), log10.(gdf.connections),
                show_empty_bins=true, color=cgrad(:plasma, scale = :exp),
                xlabel="log10 UMIs (mean: $xmean, median: $xmed)",
                ylabel="log10 connections (mean: $ymean, median: $ymed)",
                title="$R UMI Distribution", titlefont = 10, guidefont = 8,
                xlims=(0, m), ylims=(0, m), bins = num_bins)
        plot!(p2, [0, m], [0, m], color=:black, linewidth=1, legend = false)

        return p1, p2
    end
    p3, p5 = plot_umi_distributions(df, :sb1_i)
    p4, p6 = plot_umi_distributions(df, :sb2_i)

    p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
    savefig(p, joinpath(out_path, "SNR.pdf"))

    p = plot(p5, p6, layout = (2, 1), size=(7*100, 8*100))
    savefig(p, joinpath(out_path, "histograms.pdf"))

    println("done") ; flush(stdout) ; GC.gc()
    return nothing
end
