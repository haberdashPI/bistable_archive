export freqprior

struct Beta <: Stats{Bool}
  α::Float64
  β::Float64
end
update(d::Beta,x::Bool,w=1.0) = x ? Beta(d.α+w,d.β) : Beta(d.α,d.β+w)
downdate(d::Beta,x::Bool,w=1.0) = x ? Beta(d.α-w,d.β) : Beta(d.α,d.β-w)
mult(d::Beta,x) = d*x
Base.:(*)(d::Beta,x) = Beta(d.α*x,d.β*x)
Base.:(+)(a::Beta,b::Beta) = Beta(a.α+b.α,a.β+b.β)
pdf(d::Beta,x) = x ? d.α / (d.α + d.β) : d.β / (d.α + d.β)
logpdf(d::Beta,x) = log(pdf(d,x))

struct BinomialCond{T}
  at::Dict{T,Beta}
end
BinomialCond(pairs::Pair{T}...) where T = BinomialCond{T}(Dict(pairs...))
update!(d::BinomialCond,cond,x::Bool,w=1.0) =
  d.at[cond] = update(d.at[cond],x,w)
downdate!(d::BinomialCond,cond,x::Bool,w=1.0) =
  d.at[cond] = downdate(d.at[cond],x,w)
function mult!(d::BinomialCond,x)
  for cond in keys(d.at)
    d.at[cond] = mult(d.at[cond],x)
  end
end
Base.:(+)(a::BinomialCond,b::BinomialCond) =
  BinomialCond(merge(+,a.at,b.at))
function logpdf(d::BinomialCond,cond,x)
  p = d.at[cond].α / (d.at[cond].α + d.at[cond].β)
  x ? log(p) : log(1-p)
end

function freqprior(bias,N)
  Beta(N*(1-(bias/2+0.5)),N)
end

