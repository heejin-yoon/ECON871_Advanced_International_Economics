using Parameters, Random, Distributions, Plots, Printf, Interpolations, Optim, DataFrames, FixedEffectModels, LinearAlgebra, FloatingTableView, .Threads
import SpecialFunctions: erfc

rt = pwd()

include(rt * "/Problem Sets/Replicate_Ruhl-Willis(2017)/code/Ruhl-Willis(2017)_model.jl")

## Exercise b.

prim, res, sim = Initialize()
solve_model(prim, res)
simulate(prim, res)
## Check wether the model is working okay.

# data_lowestQ = DataFrame(pol_0=res.pol_func[1,:,1], pol_1=res.pol_func[2,:,1], val_0=res.val_func[1,:,1], val_1=res.val_func[2,:,1])
# writedlm("data_lowestQ.csv", Iterators.flatten(([names(data_lowestQ)], eachrow(data_lowestQ))), ',')

##

initial = [0.961, 0.047, 0.146, 0.117, 0.873]
initial = [0.961, 0.047, 0.1440402171160336, 0.12833318683668238, 0.8600259464524087]
sunkcost = Minimize_Moment(prim, res, initial)
ϕ = sunkcost.minimizer

