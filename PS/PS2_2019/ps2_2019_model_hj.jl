# ------------------------------------------------------------------------------
# Problem Set 2 (2019): The Melitz-Chaney Model
# Written by Heejin Yoon
# ------------------------------------------------------------------------------

## keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    n::Int64 = 2  # number of countries
    μ::Int64 = 1  # number of potential firms
    ρ::Float64 = 0.8  # elasticity of substitution σ = 1/(1-ρ)
    σ::Float64 = 1 / (1 - ρ)
    α::Float64 = 0.4  # consumption share of differentiated goods
    x̲::Float64 = 1.0  # lower bound of productivity
    γ::Float64 = 6.5  # Pareto distribution parameter
    L::Float64 = 1.0  # Fixed labor supply of each country
    κ::Array{Float64} = [0.1 0.2; 0.2 0.1]  # Fixed cost of to sell from i to j.
    τ::Array{Float64} = [1.0 1.15; 1.15 1.0] # Iceberg trade cost
    w::Float64 = 1.0 # wage is pinned down to be 1
end

## structure that holds model results

mutable struct Results
    Π::Array{Float64}  # profit of each country
    P::Array{Float64}  # aggregate price index of each country
    Y::Array{Float64}  # aggregate income of each country
    ϕ_cutoff::Array{Float64,2}  # cutoff productivity for each country pair
    mass_total::Array{Float64}  # mass of differentiated goods consumed
    mass_imported::Array{Float64}  # mass of differentiated imported goods consumed
    mass_domestic::Array{Float64}  # mass of differentiated domestic goods consumed
    X::Array{Float64,2}  # the value of consumption from foreign and domestic firms
end

## function for initializing model primitives and results

function Initialize()
    prim = Primitives()
    Π = ones(prim.n)
    P = ones(prim.n)
    Y = zeros(prim.n)
    ϕ_cutoff = zeros(prim.n, prim.n)
    mass_total = zeros(prim.n)
    mass_imported = zeros(prim.n)
    mass_domestic = zeros(prim.n)
    X = zeros(prim.n, prim.n)
    res = Results(Π, P, Y, ϕ_cutoff, mass_total, mass_imported, mass_domestic, X)
    prim, res
end

##

function gdp(prim::Primitives, res::Results)
    @unpack w, L, n = prim
    @unpack Π = res

    Y = zeros(n)

    for i_index = 1:n
        Y[i_index] = w * L + Π[i_index]
    end

    Y
end

## calculate the cutoff productivity of each country pair

function cutoff_productivity(prim::Primitives, res::Results)
    @unpack n, σ, w, τ, L, α, κ = prim
    @unpack Π, P, Y = res

    ϕ_cutoff = zeros(n, n)

    for i_index = 1:n
        for j_index = 1:n
            ϕ_cutoff[i_index, j_index] = (σ / (σ - 1)) * (w * τ[i_index, j_index] / P[j_index]) * (σ^(-1) * α * Y[j_index])^(-1 / (σ - 1)) * (w * κ[i_index, j_index])^(1 / (σ - 1))
        end
    end

    ϕ_cutoff
end

## calculate the price index of each country.

function price_index(prim::Primitives, res::Results)
    @unpack n, σ, γ, x̲, w, τ, μ = prim
    @unpack ϕ_cutoff = res

    P_new = zeros(n)
    f(x) = ((σ - 1) / σ)^(σ - 1) * x^(σ - γ - 2) * μ * γ * x̲^γ

    for i_index = 1:n
        for j_index = 1:n
            P_new[i_index] += (w * τ[i_index, j_index])^(1 - σ) * quadgk(f, ϕ_cutoff[i_index, j_index], Inf, rtol=1e-3)[1]
        end
        P_new[i_index] = (P_new[i_index])^(1 / (1 - σ))
    end

    P_new
end

## calculate the profit of each country.

function profit(prim::Primitives, res::Results)
    @unpack α, n, σ, γ, x̲, w, τ, μ, κ = prim
    @unpack Y, P, ϕ_cutoff = res

    Π_new = zeros(n)

    g(x) = x^(σ - γ - 2) * μ * γ * x̲^γ

    for i_index = 1:n
        for j_index = 1:n
            Π_new[i_index] += σ^(-1) * α * Y[j_index] * ((σ * w * τ[i_index, j_index]) / ((σ - 1) * P[j_index]))^(1 - σ) * quadgk(g, ϕ_cutoff[i_index, j_index], Inf, rtol=1e-3)[1] - κ[i_index, j_index] * μ * x̲^γ * ϕ_cutoff[i_index, j_index]^(-γ)
        end
    end

    Π_new
end

## define find_wage function to minimize the excess trade flow

function solve_model(prim::Primitives, res::Results)
    tol = 1e-4
    err = 1.0
    n = 1 # count
    λ = 0.5
    while true
        res.Y = gdp(prim, res)
        res.ϕ_cutoff = cutoff_productivity(prim, res)
        P_new = price_index(prim, res)
        Π_new = profit(prim, res)
        err = maximum([abs.(P_new .- res.P) abs.(Π_new .- res.Π)])
        res.P = λ .* P_new + (1 - λ) .* res.P
        res.Π = λ .* Π_new + (1 - λ) .* res.Π
        println("***** ", n, "th iteration *****")
        @printf("Absolute difference: %0.5f \n", float(err))
        println("***************************")
        n += 1
        if err < tol
            break
        end
    end
    println("Profits and price indices converged after ", n, " iterations.")
    println(" ")
    res.mass_total, res.mass_imported, res.mass_domestic = mass_goods(prim, res)
    res.X = FOB_value(prim, res)

    res
end

##

function mass_goods(prim::Primitives, res::Results)
    @unpack x̲, γ, n = prim
    @unpack ϕ_cutoff = res

    mass_imported = zeros(n)
    mass_domestic = zeros(n)

    for i_index = 1:n
        mass_imported[i_index] = 1 - (1 - x̲^γ * ϕ_cutoff[i_index, 3-i_index]^(-γ))
        mass_domestic[i_index] = 1 - (1 - x̲^γ * ϕ_cutoff[i_index, i_index]^(-γ))
    end

    mass_total = mass_imported + mass_domestic

    mass_total, mass_imported, mass_domestic
end

##

function FOB_value(prim::Primitives, res::Results)
    @unpack α, n, σ, w, τ, μ, γ, x̲ = prim
    @unpack Y, P, ϕ_cutoff = res

    X = zeros(n, n)
    h(x) = μ * γ * x̲^γ * x^(σ - γ - 2)

    for i_index = 1:n
        for j_index = 1:n
            X[i_index, j_index] = α * Y[j_index] * P[j_index]^(σ - 1) * ((σ * w * τ[i_index, j_index]) / ((σ - 1)))^(1 - σ) * quadgk(h, ϕ_cutoff[i_index, j_index], Inf, rtol=1e-3)[1]
        end
    end

    X
end
