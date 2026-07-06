using Plots
using StatsBase
using StatsPlots
using DataFrames
using KernelDensity

function filter_connections(
    df::DataFrame,
    metadata::Dict{String,Int64},
    out_path::String,
)::Nothing
    print_start("Connection filter... ")

    function connection_filter(df::DataFrame, metadata::Dict, col::Symbol, z=+3)
        R = "R"*string(col)[3] # R1 or R2

        if nrow(df) == 0
            metadata["$(R)_cxnfilter_z"] = z
            metadata["$(R)_cxnfilter_cutoff"] = 0
            metadata["$(R)_cxnfilter_beads"] = 0
            metadata["$(R)_cxnfilter"] = 0

            local p1 = plot(title="$R Connections", titlefont=10, guidefont=8,
                            xlabel="Connections (log10)", ylabel="Count")
            local p2 = plot(title="$R Max UMIs per connection", titlefont=10, guidefont=8,
                            xlabel="Maximum UMI (log10)", ylabel="Count")
            return Set{UInt64}(), p1, p2
        end

        gdf = combine(groupby(df, col), :umi => sum => :umi,
                                        :umi => length => :connections,
                                        :umi => maximum => :maximum)

        logcon = log10.(gdf.connections)
        logmax = log10.(gdf.maximum)

        # assume values above mode follow a half-normal distribution
        kdens = kde(logcon)
        xmode = kdens.x[argmax(kdens.density)]
        sd = mean(filter(x -> x >= xmode, logcon) .- xmode) * sqrt(pi) / sqrt(2)
        z3 = xmode + sd * z

        # create the label
        a = round(10^z3, digits=2)
        b = sum(gdf.connections .> 10^z3)
        c = round(b / nrow(gdf) * 100, digits=2)
        label = "z = +$z: $a\nbeads: $b ($c%)"

        # connection plot
        local p1 = barhist(logcon, bins=100, legend=false, line=0,
                           titlefont = 10, guidefont = 6, xlabelfontsize = 8, ylabelfontsize = 8,
                           title="$R Connections", xlabel="Connections (log10)", ylabel="Count")
        vline!(p1, [z3], color=:red, linestyle=:dash, linewidth=1)
        annotate!(p1, z3+Plots.xlims(p1)[2]*0.02, Plots.ylims(p1)[2]*0.95,
                  text(label, :red, 6, :left))

        # max plot
        local p2 = barhist(logmax, bins=100, legend=false, line=0,
                           titlefont = 10, guidefont = 6, xlabelfontsize = 8, ylabelfontsize = 8,
                           title="$R Max UMIs per connection", xlabel="Maximum UMI (log10)", ylabel="Count")

        gdf_remove = filter(:connections => c -> c > 10^z3, gdf)
        remove_set = Set(gdf_remove[:, col])

        metadata["$(R)_cxnfilter_z"] = z
        metadata["$(R)_cxnfilter_cutoff"] = isnan(z3) ? 0 : ceil(10^z3)
        metadata["$(R)_cxnfilter_beads"] = length(remove_set)
        metadata["$(R)_cxnfilter"] = sum(gdf_remove.umi)

        return remove_set, p1, p2
    end

    R1_remove, p1, p3 = connection_filter(df, metadata, :sb1_i)
    R2_remove, p2, p4 = connection_filter(df, metadata, :sb2_i)

    p = plot(p1, p2, p3, p4, layout = (2, 2), size=(7*100, 8*100))
    savefig(p, joinpath(out_path, "connection_filter.pdf"))

    before = sum(df.umi)
    subset!(df, :sb1_i => x -> .!in.(x, [R1_remove]), :sb2_i => y -> .!in.(y, [R2_remove]))
    after = sum(df.umi)

    metadata["umis_cxnfilter"] = before - after
    metadata["umis_final"] = sum(df.umi)
    metadata["connections_final"] = nrow(df)

    println_done()
    return nothing
end
