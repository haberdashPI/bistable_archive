using FileIO
using JLD2
using DataFrames
using Feather

include("count_lengths.jl")

dir = joinpath("..","data","count_lengths","run_2018-09-12")
# dir = joinpath(homedir(),"work","dlittle","bistable_threshold_freq","data")
all_rows = []
for_count_lengths(dir) do count_length
  push!(all_rows,DataFrame(pindex = count_length.pindex,
                           created = count_length.created,
                           ratio = count_length.ratio,
                           kind = "component"))
  push!(all_rows,DataFrame(pindex = count_length.pindex,
                           created = count_length.created,
                           ratio = count_length.bratio,
                           kind = "bandwidth"))
end
df = vcat(all_rows...)

Feather.write(joinpath("..","data","count_lengths", "run_2018-09-12",
                       "individual_results.feather"),df)

