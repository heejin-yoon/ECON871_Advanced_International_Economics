using Parameters, Plots, Printf, Distributions, DelimitedFiles, Random, Setfield, Optim, Weave

rt = pwd()

include(rt * "/Problem Sets/PS1_F2019/code/ps1_F2019_model_hj.jl")

## Exercise 1-a.

prim, res = Initialize()

T = [1.5 1.5 1.5]
res.z = draw_productivity(prim, res, T)

τ = zeros(prim.nc, prim.nc)

elapse_1a = @elapsed solve_balanced_trade(prim, res, τ)
w_1a = res.w
Π_1a = res.Π
wg_1a = welfare_gain(prim, res)

println()
println("*********** Exerise 1-a ***********")
@printf("\nIt took %0.4f seconds to solve the trade equilibrium.\n", float(elapse_1a))
@printf("\nWage of each country: [%0.4f; %0.4f; %0.4f]\n", float(w_1a[1]), float(w_1a[2]), float(w_1a[3]))
@printf("\nBilateral Trade Share: \n   Country 1: [%0.4f; %0.4f; %0.4f] \n   Country 2: [%0.4f; %0.4f; %0.4f] \n   Country 3: [%0.4f; %0.4f; %0.4f] \n\n", float(Π_1a[1, 1]), float(Π_1a[2, 1]), float(Π_1a[3, 1]), float(Π_1a[1, 2]), float(Π_1a[2, 2]), float(Π_1a[3, 2]), float(Π_1a[1, 3]), float(Π_1a[2, 3]), float(Π_1a[3, 3]))
@printf("Welfare gain of each country: [%0.4f; %0.4f; %0.4f]\n\n", float(wg_1a[1]), float(wg_1a[2]), float(wg_1a[3]))
println("***********************************\n")

## Excercise 1-b.

prim, res = Initialize()

T = [1.5 1.5 1.5]
res.z = draw_productivity(prim, res, T)

τ = zeros(prim.nc, prim.nc)
for i_index = 1:prim.nc
    for j_index = 1:prim.nc
        if i_index != j_index
            τ[i_index, j_index] = 0.1
        end
    end
end

elapse_1b = @elapsed solve_balanced_trade(prim, res, τ)
w_1b = res.w
Π_1b = res.Π
wg_1b = welfare_gain(prim, res)

println()
println("*********** Exerise 1-b ***********")
@printf("\nIt took %0.4f seconds to solve the trade equilibrium.\n", float(elapse_1b))
@printf("\nWage of each country: [%0.4f; %0.4f; %0.4f]\n", float(w_1b[1]), float(w_1b[2]), float(w_1b[3]))
@printf("\nBilateral Trade Share: \n   Country 1: [%0.4f; %0.4f; %0.4f] \n   Country 2: [%0.4f; %0.4f; %0.4f] \n   Country 3: [%0.4f; %0.4f; %0.4f] \n\n", float(Π_1b[1, 1]), float(Π_1b[2, 1]), float(Π_1b[3, 1]), float(Π_1b[1, 2]), float(Π_1b[2, 2]), float(Π_1b[3, 2]), float(Π_1b[1, 3]), float(Π_1b[2, 3]), float(Π_1b[3, 3]))
@printf("Welfare gain of each country: [%0.4f; %0.4f; %0.4f]\n\n", float(wg_1b[1]), float(wg_1b[2]), float(wg_1b[3]))
println("***********************************\n")

## Excercise 1-c.

prim, res = Initialize()

T = [1.5 1.5 1.5]
res.z = draw_productivity(prim, res, T)

τ = zeros(prim.nc, prim.nc)
τ[1, 2] = 1.05
τ[2, 1] = 1.05
τ[3, 1] = 1.3
τ[3, 2] = 1.3
τ[1, 3] = 1.3
τ[2, 3] = 1.3

elapse_1c = @elapsed solve_balanced_trade(prim, res, τ)
w_1c = res.w
Π_1c = res.Π
wg_1c = welfare_gain(prim, res)

println()
println("*********** Exerise 1-c ***********")
@printf("\nIt took %0.4f seconds to solve the trade equilibrium.\n", float(elapse_1c))
@printf("\nWage of each country: [%0.4f; %0.4f; %0.4f]\n", float(w_1c[1]), float(w_1c[2]), float(w_1c[3]))
@printf("\nBilateral Trade Share: \n   Country 1: [%0.4f; %0.4f; %0.4f] \n   Country 2: [%0.4f; %0.4f; %0.4f] \n   Country 3: [%0.4f; %0.4f; %0.4f] \n\n", float(Π_1c[1, 1]), float(Π_1c[2, 1]), float(Π_1c[3, 1]), float(Π_1c[1, 2]), float(Π_1c[2, 2]), float(Π_1c[3, 2]), float(Π_1c[1, 3]), float(Π_1c[2, 3]), float(Π_1c[3, 3]))
@printf("Welfare gain of each country: [%0.4f; %0.4f; %0.4f]\n\n", float(wg_1c[1]), float(wg_1c[2]), float(wg_1c[3]))
println("***********************************\n")

## Excercise 1-d.

prim, res = Initialize()

T = [1.5 3.0 1.5]
res.z = draw_productivity(prim, res, T)

τ = zeros(prim.nc, prim.nc)
τ[1, 2] = 1.05
τ[2, 1] = 1.05
τ[3, 1] = 1.3
τ[3, 2] = 1.3
τ[1, 3] = 1.3
τ[2, 3] = 1.3

elapse_1d = @elapsed solve_balanced_trade(prim, res, τ)
w_1d = res.w
Π_1d = res.Π
wg_1d = welfare_gain(prim, res)

println()
println("*********** Exerise 1-d ***********")
@printf("\nIt took %0.4f seconds to solve the trade equilibrium.\n", float(elapse_1d))
@printf("\nWage of each country: [%0.4f; %0.4f; %0.4f]\n", float(w_1d[1]), float(w_1d[2]), float(w_1d[3]))
@printf("\nBilateral Trade Share: \n   Country 1: [%0.4f; %0.4f; %0.4f] \n   Country 2: [%0.4f; %0.4f; %0.4f] \n   Country 3: [%0.4f; %0.4f; %0.4f] \n\n", float(Π_1d[1, 1]), float(Π_1d[2, 1]), float(Π_1d[3, 1]), float(Π_1d[1, 2]), float(Π_1d[2, 2]), float(Π_1d[3, 2]), float(Π_1d[1, 3]), float(Π_1d[2, 3]), float(Π_1d[3, 3]))
@printf("Welfare gain of each country: [%0.4f; %0.4f; %0.4f]\n\n", float(wg_1d[1]), float(wg_1d[2]), float(wg_1d[3]))
println("***********************************\n")

## Excercise 1-e.

prim, res = Initialize()

prim = @set prim.θ = 8.0
T = [1.5 1.5 1.5]
res.z = draw_productivity(prim, res, T)

τ = zeros(prim.nc, prim.nc)
τ[1, 2] = 1.05
τ[2, 1] = 1.05
τ[3, 1] = 1.3
τ[3, 2] = 1.3
τ[1, 3] = 1.3
τ[2, 3] = 1.3

elapse_1e = @elapsed solve_balanced_trade(prim, res, τ)
w_1e = res.w
Π_1e = res.Π
wg_1e = welfare_gain(prim, res)

println()
println("*********** Exerise 1-e ***********")
@printf("\nIt took %0.4f seconds to solve the trade equilibrium.\n", float(elapse_1e))
@printf("\nWage of each country: [%0.4f; %0.4f; %0.4f]\n", float(w_1e[1]), float(w_1e[2]), float(w_1e[3]))
@printf("\nBilateral Trade Share: \n   Country 1: [%0.4f; %0.4f; %0.4f] \n   Country 2: [%0.4f; %0.4f; %0.4f] \n   Country 3: [%0.4f; %0.4f; %0.4f] \n\n", float(Π_1e[1, 1]), float(Π_1e[2, 1]), float(Π_1e[3, 1]), float(Π_1e[1, 2]), float(Π_1e[2, 2]), float(Π_1e[3, 2]), float(Π_1e[1, 3]), float(Π_1e[2, 3]), float(Π_1e[3, 3]))
@printf("Welfare gain of each country: [%0.4f; %0.4f; %0.4f]\n\n", float(wg_1e[1]), float(wg_1e[2]), float(wg_1e[3]))
println("***********************************\n")

println("\nAll done!")

##

weave("ps1_2019_model_hj.jl", doctype="md2pdf", out_path="weave")
