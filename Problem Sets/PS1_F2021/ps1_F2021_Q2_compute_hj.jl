using Parameters, Random, Optim, Distributions, Plots, Printf, Interpolations, DataFrames, FixedEffectModels, Weave

rt = pwd()

include(rt * "/Problem Sets/PS1_F2021/ps1_F2021_Q2_model_hj.jl")

## Q2. The sunk-cost model

# (a)-1. solve for the stationary equilibrium 

prim, res, sim = Initialize()
solve_model(prim, res)

# (a)-2. simulate to calculate the key moment

sim.moments_sim = simulation(prim, res, sim)

println("*********** Exerise 2-(a) ***********")
@printf("\ni. Export participation rate: %0.4f\n", float(sim.moments_sim[1]))
@printf("\nii. Export starter rate: %0.4f\n", float(sim.moments_sim[2]))
@printf("\niii. Export stopper rate: %0.4f\n", float(sim.moments_sim[3]))
@printf("\niv.Mean export intensity of exporters: %0.4f\n", float(sim.moments_sim[4]))
@printf("\nv. Autocorrelation of domestic sales: %0.4f\n", float(sim.moments_sim[5]))
println("***********************************\n")


## (b). plot the polcy function

plot(res.a_grid, [res.pol_func[2, :] res.pol_func[1, :]], xlims=(0, 3), ylims=(0, 1.2), yticks=(0:1:1), title="Policy Functions", xlabel="a values", ylabel="X prime", labels=["X = 1" "X = 0"], legend=:topleft)

##

println("\nAll done!")

weave("ps1_F2021_Q2_model_hj.jl", doctype = "md2pdf", out_path = "weave")
