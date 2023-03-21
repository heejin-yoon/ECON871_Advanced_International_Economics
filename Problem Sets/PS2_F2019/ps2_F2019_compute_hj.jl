using Parameters, Plots, Printf, DelimitedFiles, Setfield, QuadGK, Weave

rt = pwd()

include(rt * "/Problem Sets/PS2_F2019/ps2_F2019_model_hj.jl")

## Exercise b.

prim, res = Initialize()
results_b = solve_model(prim, res)

## Exercise c.

prim, res = Initialize()
τ_c = [1.0 1.0; 1.0 1.0]
prim = @set prim.τ = τ_c
results_c = solve_model(prim, res)

## Exercise d.

elasticity = zeros(3)

prim, res = Initialize()
τ_b = [1.0 1.15; 1.15 1.0]
τ_d1 = [1.0 1.14999; 1.14999 1.0]
prim = @set prim.τ = τ_d1
results_d1 = solve_model(prim, res)
elasticity[1] = -((results_d1.X[2, 1] / results_d1.X[1, 1]) / (results_b.X[2, 1] / results_b.X[1, 1]) - 1) / (τ_d1[2, 1] / τ_b[2, 1] - 1)

prim, res = Initialize()
τ_d2 = [1.0 1.05; 1.05 1.0]
prim = @set prim.τ = τ_d2
results_d2 = solve_model(prim, res)
elasticity[2] = -((results_d2.X[2, 1] / results_d2.X[1, 1]) / (results_b.X[2, 1] / results_b.X[1, 1]) - 1) / (τ_d2[2, 1] / τ_b[2, 1] - 1)

prim, res = Initialize()
τ_d3 = [1.0 1.0; 1.0 1.0]
prim = @set prim.τ = τ_d3
results_d3 = solve_model(prim, res)
elasticity[3] = -((results_d3.X[2, 1] / results_d3.X[1, 1]) / (results_b.X[2, 1] / results_b.X[1, 1]) - 1) / (τ_d3[2, 1] / τ_b[2, 1] - 1)

##

println("\nAll done!")

weave("ps2_2019_model_hj.jl", doctype="md2pdf", out_path="weave")
