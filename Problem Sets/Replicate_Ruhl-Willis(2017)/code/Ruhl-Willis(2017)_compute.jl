using Parameters, Random, Distributions, Plots, Printf, Interpolations, Optim, DataFrames, FixedEffectModels, LinearAlgebra, FloatingTableView, .Threads
import SpecialFunctions: erfc

rt = pwd()

include(rt * "/Problem Sets/Replicate_Ruhl-Willis(2017)/code/Ruhl-Willis(2017)_model.jl")

prim, res, sim = Initialize()


## Check whether the model is working okay

solve_model(prim, res)
sim.moments_sim = simulate(prim, res)
diff = abs.(sim.moments_sim - prim.moments_data) 
moments = DataFrame([["Starter Rate", "Stopper Rate", "Average Export-Sales Ratio", "Coefficient of Variation", "Slope of Domestic Sales Regression"] prim.moments_data sim.moments_sim diff], [:Moment, :Data, :Simulation_Baseline, :Difference])
show(stdout, moments)


## Sunk-Cost Model (Baseline)

initial = [0.961, 0.047, 0.146, 0.117, 0.873]
sunkcost = Minimize_Moment(prim, res, initial)
ϕ_sunkcost = sunkcost.minimizer

