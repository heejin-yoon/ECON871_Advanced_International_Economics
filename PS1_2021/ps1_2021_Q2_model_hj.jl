# -------------------------------------------------------------------------------
# Problem Set 1 (2021): Q2. The Sunk Cost Model
# Written by Heejin Yoon
# -------------------------------------------------------------------------------

## keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    θ::Float64 = 4.00     # elasticity of substitution
    ā::Float64 = 1.00     # mean of log(a) process
    β::Float64 = 0.90     # discount rate
    τ::Float64 = 1.10     # tariff
    w::Float64 = 1.14     # (patial) equilibrium wage
    ρ::Float64 = 0.92     # log(a) process parameter
    σ::Float64 = 0.25     # log(a) process parameter
    ξ::Float64 = 1.20     # variable export cost
    f₀::Float64 = 0.40    # entry cost
    f₁::Float64 = 0.28    # continuing cost
    n_a::Int64 = 1000     # number of a_grid
    T::Int64 = 1000       # number of periods for simulation
    N::Int64 = 15000      # number of firms for simulation
    drop::Int64 = 200     # periods that will be dropped
    # moments_data::Array{Float64} = [] # needed if we want to estimate key parameters using SMM.
end

## structure that holds model results

mutable struct Results
    Π_a::Array{Float64,2}
    a_grid::Array{Float64}
    val_func::Array{Float64,2}
    pol_func::Array{Int64,2}
end

## structure that holds simulation outcomes

mutable struct Simulations
    a_sim::Array{Float64,2}
    pol_func_sim::Array{Float64,2}
    R_domestic_sim::Array{Float64,2}
    R_exports_sim::Array{Float64,2}
    moments_sim::Array{Float64}
end

## function for initializing model primitives and results

function Initialize()
    prim = Primitives()
    Π_a = tauchen_method(prim.n_a, prim.ρ, prim.σ, log(prim.ā))[1]
    a_grid = exp.(tauchen_method(prim.n_a, prim.ρ, prim.σ, log(prim.ā))[2])
    val_func = zeros(2, prim.n_a)
    pol_func = zeros(Int64, 2, prim.n_a)
    a_sim = zeros(prim.T - prim.drop, prim.N)
    pol_func_sim = zeros(Int64, prim.T - prim.drop, prim.N)
    R_domestic_sim = zeros(prim.T - prim.drop, prim.N)
    R_exports_sim = zeros(prim.T - prim.drop, prim.N)
    moments_sim = zeros(5)
    res = Results(Π_a, a_grid, val_func, pol_func)
    sim = Simulations(a_sim, pol_func_sim, R_domestic_sim, R_exports_sim, moments_sim)
    prim, res, sim
end

## discretize AR(1) process using Tauchen's method (instead of manually constructing it, we can use QuantEcon.tauchen() function)

function tauchen_method(N::Int64, ρ::Float64, σ::Float64, μ::Float64)

    dist = Normal(0, σ)
    σ_y = sqrt(σ^2 / (1 - ρ^2))

    y_bar = 3 * σ_y # grid ranges 3 * standard deviation of y
    y = collect(range(-y_bar, stop=y_bar, length=N))
    s = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(N, N)
    for row_index = 1:N
        Π[row_index, 1] = cdf(dist, y[1] - ρ * y[row_index] + s / 2)
        Π[row_index, N] = 1 - cdf(dist, y[N] - ρ * y[row_index] + s / 2)
        for col_index = 2:N-1
            Π[row_index, col_index] = cdf(dist, y[col_index] - ρ * y[row_index] + s / 2) - cdf(dist, y[col_index] - ρ * y[row_index] - s / 2)
        end
    end

    Π = Π ./ sum(Π, dims=2)
    y = y .+ μ

    Π, y
end

## fixed cost depending on X and X_prime values.

function fixed_cost(prim::Primitives, res::Results, X::Int64, X_prime::Int64)
    @unpack f₀, f₁ = prim

    if X_prime == 0
        f = 0.0
    elseif X_prime == 1
        if X == 0
            f = f₀
        elseif X == 1
            f = f₁
        end
    end

    f
end

## calculate the profit function.

function profit_function(prim::Primitives, res::Results, X::Int64, a_value::Float64)
    @unpack θ, ξ, w, τ = prim

    α = (θ - 1) / θ
    # K = (a_value/((θ-1)*(ξ-1)))^((θ-1)/θ) * (1 + X * (θ*ξ - θ - ξ)/ξ)
    # K = ((a_value/((θ-1)*(ξ-1)))^((θ-1)/θ)*(1-X/ξ) + X/ξ * (a_value/((θ-1)*(ξ-1)))^(-1/θ)*a_value)
    # K = a_value/(1 + ξ * τ^(-θ)) * (1 + X * τ^(-θ))
    K = ((τ^θ / ξ) / ((τ^θ / ξ) + 1) * a_value)^α + X / τ * ((a_value / ξ) / ((τ^θ / ξ) + 1))^α
    Π = (1 - α) * K^(1 / (1 - α)) * (α / w)^(α / (1 - α))

    Π
end

## calculate the profit of each country.

function Bellman(prim::Primitives, res::Results)
    @unpack β, n_a = prim
    @unpack val_func, pol_func, Π_a, a_grid = res

    val_func_new = zeros(2, n_a)
    pol_func_new = zeros(Int64, 2, n_a)

    for a_index = 1:n_a
        Π_0 = profit_function(prim, res, 0, a_grid[a_index])
        Π_1 = profit_function(prim, res, 1, a_grid[a_index])

        for X_index = 1:2
            val_func_new_ex = Π_1 - fixed_cost(prim, res, X_index - 1, 1)
            val_func_new_no_ex = Π_0
            for ap_index = 1:n_a
                val_func_new_ex += β * Π_a[a_index, ap_index] * val_func[2, ap_index]
                val_func_new_no_ex += β * Π_a[a_index, ap_index] * val_func[1, ap_index]
            end

            if val_func_new_ex > val_func_new_no_ex
                pol_func_new[X_index, a_index] = 1
                val_func_new[X_index, a_index] = val_func_new_ex
            else
                pol_func_new[X_index, a_index] = 0
                val_func_new[X_index, a_index] = val_func_new_no_ex
            end
        end
    end

    val_func_new, pol_func_new
end

## define find_wage function to minimize the excess trade flow

function solve_model(prim::Primitives, res::Results)
    tol = 1e-4
    err = 1.0
    n = 1

    while true
        val_func_new, pol_func_new = Bellman(prim, res)
        err = maximum(abs.(res.val_func - val_func_new))
        res.val_func = val_func_new
        res.pol_func = pol_func_new
        println("***** ", n, "th iteration *****")
        @printf("Absolute difference: %0.5f \n", float(err))
        println("***************************")
        n += 1
        if err < tol
            break
        end
    end
    println("Value function converged after ", n, " iterations.")
    println(" ")

end

## get index from a grid

function get_index(val::Float64, grid::Array{Float64,1})

    interp = interpolate(grid, BSpline(Linear()))
    find_index(k) = abs(interp(k) - val)
    index = optimize(find_index, 1.0, length(grid)).minimizer
    index = round(Int, index)

    index
end

## simulate productivty shock

function simulate_a(prim::Primitives, res::Results)
    @unpack ā, ρ, σ, T, N, drop = prim

    Random.seed!(12341234)
    log_a_sim = zeros(T, N)

    for N_index = 1:N
        η = randn() * sqrt(σ^2 / (1 - ρ^2))
        log_a_sim[1, N_index] = η
        for T_index = 2:T
            log_a_sim[T_index, N_index] = (1 - ρ) * log(ā) + ρ * log_a_sim[T_index-1, N_index] + randn() * σ
        end
    end
    a_sim = exp.(log_a_sim)

    a_sim
end

## calculate the revenue function: something is wrong.

function revenue_function(prim::Primitives, res::Results, X::Int64, a_value::Float64)
    @unpack θ, ξ, w, τ = prim

    α = (θ - 1) / θ
    K = ((τ^θ / ξ) / ((τ^θ / ξ) + 1) * a_value)^α + X / τ * ((a_value / ξ) / ((τ^θ / ξ) + 1))^α

    Y = K^(α / (1 - α)) * (α / w)^(α / (1 - α))

    R_domestic = ((τ^θ / ξ) / ((τ^θ / ξ) + 1) * a_value)^α * Y
    R_exports = X / τ * ((a_value / ξ) / ((τ^θ / ξ) + 1))^α * Y

    R_domestic, R_exports

end


## get policy function, domestic and export sales revenues using simulated productivity a

function simulate_pol_rev(prim::Primitives, res::Results, sim::Simulations)
    @unpack N, T = prim
    @unpack pol_func, a_grid = res
    @unpack a_sim = sim

    pol_func_sim = zeros(Int64, T, N)
    R_domestic_sim = zeros(T, N)
    R_exports_sim = zeros(T, N)

    for N_index = 1:N
        for T_index = 1:T
            a_index = get_index(a_sim[T_index, N_index], a_grid)
            a_value = a_grid[a_index]
            if T_index == 1
                pol_func_sim[T_index, N_index] = 0
            else
                pol_func_sim[T_index, N_index] = pol_func[pol_func_sim[T_index-1, N_index]+1, a_index]
            end
            R_domestic_sim[T_index, N_index], R_exports_sim[T_index, N_index] = revenue_function(prim, res, pol_func_sim[T_index, N_index], a_value)
        end
    end

    pol_func_sim, R_domestic_sim, R_exports_sim
end

##

function simulation(prim::Primitives, res::Results, sim::Simulations)
    @unpack T, N, drop = prim

    sim.a_sim = simulate_a(prim, res)
    sim.pol_func_sim, sim.R_domestic_sim, sim.R_exports_sim = simulate_pol_rev(prim, res, sim)

    X = sim.pol_func_sim[1, :]
    pol_func_sim = sim.pol_func_sim
    export_part_rate = mean(pol_func_sim[drop+1:T, :])

    starterrate = zeros(T)
    stopperrate = zeros(T)

    for T_index = 2:T
        total_nonexperter = N - sum(X)
        X_prime = pol_func_sim[T_index, :]
        starter = 0
        stopper = 0

        for N_index = 1:N
            if X[N_index] == 0
                if X_prime[N_index] == 1
                    starter = starter + 1
                end
            else
                if X_prime[N_index] == 0
                    stopper = stopper + 1
                end
            end
        end
        starterrate[T_index] = starter / total_nonexperter
        stopperrate[T_index] = stopper / (N - total_nonexperter)
        X = X_prime
    end

    avg_starterrate = mean(starterrate[drop+1:T])
    avg_stopperrate = mean(stopperrate[drop+1:T])
    plot([starterrate[drop+1:T], avg_starterrate * ones(T - drop)], labels=["Starter Rate" "LR Average"], legend=:topright)
    plot([stopperrate[drop+1:T], avg_stopperrate * ones(T - drop)], labels=["Stopper Rate" "LR Average"], legend=:topright)

    R_domestic_sim = sim.R_domestic_sim[drop+1:T, :]
    R_exports_sim = sim.R_exports_sim[drop+1:T, :]
    Y = T - drop

    exportsales = zeros(Y, N)

    for N_index = 1:N
        for Y_index = 1:Y
            exportsales[Y_index, N_index] = R_exports_sim[Y_index, N_index] / (R_domestic_sim[Y_index, N_index] + R_exports_sim[Y_index, N_index])
        end
    end

    exporter = (R_exports_sim .> 0)
    export_intensity = sum(exportsales) / (sum(exporter))

    y = zeros(N * Y)
    y_lag = zeros(N * Y)
    firmno = zeros(N * Y)
    year = zeros(N * Y)

    for N_index = 1:N
        for Y_index = 1:Y
            y[(N_index-1)*Y+Y_index] = R_domestic_sim[Y_index, N_index]
            firmno[(N_index-1)*Y+Y_index] = N_index
            year[(N_index-1)*Y+Y_index] = Y_index
        end
        for Y_index = 2:Y
            y_lag[(N_index-1)*Y+Y_index] = R_domestic_sim[Y_index-1, N_index]
        end
        y_lag[(N_index-1)*Y+1] = NaN
    end

    df = DataFrame(y_lag=y_lag, y=y, firmno=firmno, year=year)
    df = filter(row -> !isnan(row.y_lag), df)

    γ = reg(df, @formula(y ~ y_lag + fe(firmno) + fe(year))).coef[1]

    [export_part_rate, avg_starterrate, avg_stopperrate, export_intensity, γ]
end

##

# function m_hat(f₀, f₁, C_star, σ_ϵ, ρ_ϵ)
#     res.f₀ = f₀
#     res.f₁ = f₁
#     res.C_star = C_star
#     res.σ_ϵ = σ_ϵ
#     res.ρ_ϵ = ρ_ϵ

#     solve_model(prim, res)
#     simulate(prim, res)
# end

# ##

# function MoM(prim::Primitives, res::Results, sim::Simulations)
#     @unpack f₀, f₁, C_star, σ_ϵ, ρ_ϵ = res
#     @unpack moments_sim = sim

#     W = Matrix{Float64}(I, 5, 5)
#     J(f₀, f₁, C_star, σ_ϵ, ρ_ϵ) = (prim.moments_data - m_hat(f₀, f₁, C_star, σ_ϵ, ρ_ϵ))' * W * (prim.moments_data - m_hat(f₀, f₁, C_star, σ_ϵ, ρ_ϵ))
#     lower = [0.1, 0.0, 0.0, 0.0, 0.3]
#     upper = [1.0, 0.5, 0.6, 0.3, 1.0]
#     # b_initial = (lower + upper)./2
#     b_hat = optimize(b -> J(b[1], b[2], b[3], b[4], b[5]), lower, upper, [f₀, f₁, C_star, σ_ϵ, ρ_ϵ]).minimizer

#     b_hat
# end

##
