using CodecZlib

function write_barcode_set(sb_is::Set{UInt64}, path::String)::Nothing
    open(GzipCompressorStream, path, "w") do io
        print(io, "") # ensure that an empty gz file is created if the tab is empty
        for sb_i in sb_is
            println(io, sb_i)
        end
    end
    return nothing
end

function read_barcode_set(path::String)::Set{UInt64}
    sb_is = Set{UInt64}()
    open(GzipDecompressorStream, path, "r") do io
        for line in eachline(io)
            push!(sb_is, parse(UInt64, line))
        end
    end
    return sb_is
end
