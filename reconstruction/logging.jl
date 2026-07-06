using Dates

# Used to report the elapsed time plus current memory usage when printing "done" messages.
# Unlike the Profile module, the memory reported is not the total memory allocated during the execution,
# but rather the current live memory usage at the time of printing after performing garbage collection.
# For actual profiling see https://docs.julialang.org/en/v1/manual/profile/

const last_print = Ref{DateTime}(now())
const CONCISE_ABBREVS = Dict{DataType,String}(
    Year => "y",
    Month => "mo",
    Week => "w",
    Day => "d",
    Hour => "h",
    Minute => "m",
    Second => "s",
    Millisecond => "ms",
    Microsecond => "µs",
    Nanosecond => "ns",
)

function concise_show(period::Period)::String
    secs = Dates.canonicalize(ceil(period, Second))
    if isempty(secs.periods)
        return "0 s"
    end
    parts = ["$(p.value) $(CONCISE_ABBREVS[typeof(p)])" for p in secs.periods]
    return join(parts, ", ")
end

function print_start(msg::String)::Nothing
    print(msg)
    flush(stdout)
    last_print[] = now()
    return nothing
end

function println_done(msg::String = "done")::Nothing
    elapsed = now() - last_print[]
    GC.gc()
    mem_gb = ceil(Int64, Base.gc_live_bytes() / 1024^3)
    println("$msg ($(concise_show(elapsed)), $mem_gb GB)")
    flush(stdout)
    last_print[] = now()
    return nothing
end
