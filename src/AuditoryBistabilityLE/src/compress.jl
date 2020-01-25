using CodecZlib
using Statistics
export compress, decompress

struct CompressedMask{A}
  size::Tuple{Int,Int}
  axes::A
  data::Array{UInt8}
end

# NOTE: this somewhat odd duplication of types is to deal with legacy code. I
# want the newly stored data to use only builtin types of julia, which save to
# jld2 files in a simple format. This helps ensure it can be readily loaded in
# with future versions of julia and in other languages

struct CMask{T,F}
  size::Tuple{Int,Int}
  data::Array{UInt8}
  times_s::T
  freqs_Hz::F
end

function CMask(x::CompressedMask)
  @assert x.axes[1] isa Axis{:time}
  @assert x.axes[2] isa Axis{:freq}
  @assert unit(eltype(x.axes[1])) == s
  @assert unit(eltype(x.axes[2])) == Hz

  CMask(x.size,x.data,ustrip.(x.axes[1].val),
                         ustrip.(x.axes[2].val))
end

function CompressedMask(x::CMask)
  CompressedMask(x.size, (Axis{:time}(s .* x.times_s),
                          Axis{:freq}(Hz .* x.freqs_Hz)), x.data)
end

compress_axes(x::MetaUnion{AxisArray}) = AxisArrays.axes(x)
compress_axes(x) = nothing

function compress(x::AbstractMatrix)
  quantized = Array{UInt8}(undef,size(x))
  @. quantized = floor(UInt8,max(0,x/$maximum(x))*typemax(UInt8))
  CMask(CompressedMask(size(x), compress_axes(x),
                                        transcode(ZlibCompressor,vec(quantized))))
end

withaxes(x,axes) = AxisArray(x, axes...)
withaxes(x,::Nothing) = x

decompress(x::CMask) = decompress(CompressedMask(x))
function decompress(x::CompressedMask)
  mask = transcode(ZlibDecompressor,x.data) ./ typemax(UInt8)
  withaxes(reshape(mask,x.size...), x.axes)
end

function CorticalSpectralTemporalResponses.audiospect(x::CompressedMask,settings)
  settings = read_settings(settings)
  audiospect(decompress(x);settings.freqs.analyze...)
end

