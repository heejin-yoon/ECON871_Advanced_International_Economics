using Parameters, Random, Distributions, Plots, Printf, Interpolations, Optim, DataFrames, FixedEffectModels, LinearAlgebra
import SpecialFunctions: erfc

cd("C:/Users/hyoon76/OneDrive - UW-Madison/3.Wisconsin/2022 Fall/Econ871/Problem Sets/Ruhl-Willis(2017)/")
include("Ruhl-Willis(2017)_hj.jl")

## Exercise b.

prim, res, sim = Initialize()
solve_model(prim, res)

## Check wether the model is working okay.

data_lowestQ = DataFrame(pol_0=res.pol_func[1,:,1], pol_1=res.pol_func[2,:,1], val_0=res.val_func[1,:,1], val_1=res.val_func[2,:,1])
writedlm("data_lowestQ.csv", Iterators.flatten(([names(data_lowestQ)], eachrow(data_lowestQ))), ',')

##

sim.moments_sim = simulate(prim, res)
