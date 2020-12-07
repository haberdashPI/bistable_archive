
totime(x) = x.*ms
tofreq(x) = x.*Hz
const data_dir = if occursin("Mycroft",gethostname())
  joinpath(@__DIR__,"..","data")
else
  joinpath(homedir(),"work","dlittle","bistable_individual")
end

from_unit_table = Dict(r"τ" => totime, r"Δt" => totime, r"^f$" => tofreq)
function handle_units!(df)
  for col in names(df)
    for (pattern,fn) in pairs(from_unit_table)
      if occursin(pattern,string(col))
        df[col] = fn.(df[col])
      end
    end
  end
  df
end

function getparams(filterfn,file)
  asdict(params) = Dict(k => params[1,k] for k in names(params))

  params = load_params(file)
  rows = filter(ri -> filterfn(ri,asdict(params[ri,:])),Base.axes(params,1))
  if length(rows) > 1
    @warn("Specification ambiguous, multiple rows match.")
  elseif length(rows) < 1
    error("No rows matching specification")
  else
    asdict(params[rows[1],:])
  end
end

function withunit(x,unitname)
  if unitname == "nothing"
    x
  else
    x * @eval(@u_str($unitname))
  end
end

# handle all of the legacy formats I have used to store
# data
function load_params(params)
  # version 0.1.0: save as a feather file
  # (move away form this because feather files don't currently guarantee
  # that new versions will be compatible with the format)
  if occursin(r"\.feather$",params)
    @warn "Using old parameter file format, run `clean_up_param_files.jl`"
    handle_units!(Feather.read(params))

  elseif occursin(r"\.jld2$",params)
    data = jldopen(params,"r")
    # version 0.2.0: save dataframe to jlds
    # (move away form this because DataFrames has no gaurantee of consistency
    # across versions yet)
    if "version" ∉ keys(data) # version 0.2.0
      @warn "Using old parameter file format, run `clean_up_param_files.jl`"
      data["params"]

    # version 0.3.0: save as arrays of columns with any units stored separately
    # this way, just basic types are stored, and hdf is a well established
    # standard
    elseif v"0.3" <= data["version"] < v"0.4"
      df = DataFrame()
      for col in keys(data["params"])
        df[:,Symbol(col)] = withunit(data["params"][col],data["units"][col])
      end
      df
    end
  end
end

units = [:ms,:s,:Hz,:kHz,:cycoct]
for unit in units
  @eval unitname(x::Quantity{<:Any,<:Any,typeof($unit)}) = $(string(unit))
end
unitname(x::Number) = "nothing"

const param_file_version = v"0.3.0"
function save_params(file,params)
  jldopen(file,"w") do data
    data["version"] = param_file_version
    for (name,val) in eachcol(params)
      data["params/"*string(name)] = ustrip.(val)
      data["units/"*string(name)] = unitname(first(val))
    end
  end
end

