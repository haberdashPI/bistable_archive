include(joinpath(@__DIR__,"interactive_setup.jl"))
curplotdir = joinpath(plotdir,string(Date(now())))
isdir(curplotdir) || mkdir(curplotdir)

datadir = joinpath(@__DIR__,"..","data","buildup","2019-08-01")
files = readdir(datadir)

df = mapreduce(vcat,files) do file
  if !isfile(joinpath(datadir,file))
    return DataFrame()
  end
  df = DataFrame(CSV.File(joinpath(datadir,file)))
  m = match(r"^([a-z]+)_df([0-9]+)_([0-9-]+)\.csv",file)
  if isnothing(m)
    @warn("Filename $file doesn't match the expected naming convention.")
  else
    df[!,:level] .= titlecase(m[1])
    df[!,:df] .= parse(Int,m[2])
    df[!,:date] .= Date(m[3])
    df
  end
end

means = by(df,[:level,:df]) do df
  buildup_mean(df,delta=0.3,length=12)
end

R"""
library(dplyr)
df = $means
df$level = factor(df$level, levels=c("Peripheral","Central","Object","Combined"), ordered=T)
df = df %>% arrange(level)
ggplot(df,aes(x=time,y=value,color=factor(df))) + geom_line() +
  facet_grid(~level) + xlim(0,10) + theme_classic()
ggsave($(joinpath(curplotdir,"buildup.pdf")))
"""
