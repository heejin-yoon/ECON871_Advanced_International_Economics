## keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    r::Float64 = 0.109  # annual interest rate
    R::Float64 = (1 + r)^(1/4)-1  # quarterly interest rate
    ρ_Q::Float64 = 0.826  # real effective exchange rate process
    σ_Q::Float64 = 0.036  # real effective exchange rate process
    n_Q::Int64 = 11
    Π_Q::Array{Float64, 2} = tauchen(n_Q, ρ_Q, σ_Q, 3.0)[1]
    Q_grid::Array{Float64} = tauchen(n_Q, ρ_Q, σ_Q, 3.0)[2]
    α_N::Float64 = 0.45  # labor share
    α_K::Float64 = 0.55  # capital share
    l::Int64 = 63  # median plant employee
    θ::Float64 = 5.0  # elasticity of substitution
    T::Int64 = 420
    N::Int64 = 1914
    drop::Int64 = 20
    moments_data::Array{Float64} = [0.0517; 0.1062; 0.1346; 0.2090; 0.6482]
end

## structure that holds model results

mutable struct Results
    f₀::Float64  # entry cost
    f₁::Float64  # continuing cost
    C_star::Float64  # world aggregate demand
    σ_ϵ::Float64  # ϵ process
    ρ_ϵ::Float64  # ϵ process
    n_ϵ::Int64
    Π_ϵ::Array{Float64, 2}
    ϵ_grid::Array{Float64}
    val_func::Array{Float64, 3}
    pol_func::Array{Int64, 3}
end

## structure that holds model results

mutable struct Simulations
    Q_sim::Array{Float64}
    ϵ_sim::Array{Float64, 2}
    pol_func_sim::Array{Float64, 2}
    domestic_sim::Array{Float64, 2}
    exports_sim::Array{Float64, 2}
    moments_sim::Array{Float64}
end

## function for initializing model primitives and results

function Initialize()
    prim = Primitives()
    f₀ = 0.961  # entry cost
    f₁ = 0.047  # continuing cost
    C_star = 0.146  # world aggregate demand
    σ_ϵ = 0.116  # ϵ process
    ρ_ϵ = 0.873  # ϵ process
    n_ϵ = 455
    Π_ϵ = tauchen(n_ϵ, ρ_ϵ, σ_ϵ, 3.0)[1]
    ϵ_grid = tauchen(n_ϵ, ρ_ϵ, σ_ϵ, 3.0)[2]
    val_func = zeros(2, n_ϵ, prim.n_Q)
    pol_func = zeros(Int64, 2, n_ϵ, prim.n_Q)
    Q_sim = zeros(prim.T - prim.drop)
    ϵ_sim = zeros(prim.T - prim.drop, prim.N)
    pol_func_sim = zeros(prim.T - prim.drop, prim.N)
    domestic_sim = zeros(prim.T - prim.drop, prim.N)
    exports_sim = zeros(prim.T - prim.drop, prim.N)
    moments_sim = zeros(length(prim.moments_data))
    res = Results(f₀, f₁, C_star, σ_ϵ, ρ_ϵ, n_ϵ, Π_ϵ, ϵ_grid, val_func, pol_func)
    sim = Simulations(Q_sim, ϵ_sim, pol_func_sim, domestic_sim, exports_sim, moments_sim)
    prim, res, sim
end

##

function tauchen(N::Int64, ρ::Float64, σ::Float64, μ::Float64)
    # Get discretized space using Tauchen's method (https://julia.quantecon.org/tools_and_techniques/finite_markov.html)

    dist = Normal(0, σ)
    σ_y = sqrt(σ^2 / (1 - ρ^2))

    y_bar = 3*σ_y
    y = collect(range(-y_bar, stop = y_bar, length = N))
    s = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(N, N)
    for row_index = 1:N
        Π[row_index, 1] = cdf(dist, y[1] - ρ*y[row_index] + s/2)
        Π[row_index, N] = 1 - cdf(dist, y[N] - ρ*y[row_index] + s/2)
        for col_index = 2:N-1
            Π[row_index, col_index] = cdf(dist, y[col_index] - ρ*y[row_index] + s/2) - cdf(dist, y[col_index] - ρ*y[row_index] - s/2)
        end
    end

    Π = Π./sum(Π, dims = 2)
    y = y .+ μ
    Π, y
end



## calculate the cutoff productivity of each country pair

function fixed_cost(prim::Primitives, res::Results, X::Int64, X_prime::Int64)
    @unpack f₀, f₁ = res

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

##

function plant_revenue(prim::Primitives, res::Results, X_prime::Int64, Q_value::Float64, ϵ_value::Float64)
    @unpack θ, α_N, α_K, l, r = prim
    @unpack C_star = res

    w = α_N *((θ-1)/θ) * l^(α_N *(θ-1)/θ-1)
    X = (1 + X_prime * (exp(Q_value)^θ) * C_star)^(1/θ) * exp(ϵ_value)

    a = α_N * (θ-1)/θ
    b = α_K * (θ - 1)/θ
    Y = exp(ϵ_value) * X^((a + b)/(1-a-b)) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b));
    domestic = (1 / (X_prime * exp(Q_value)^θ * C_star + 1))^((θ - 1)/θ) * Y;
    exports = X_prime * exp(Q_value) * C_star^(1/θ) * (1/(1 + exp(Q_value)^(-θ)/C_star))^((θ - 1)/θ) * Y;

    Y, domestic, exports
end

## calculate the price index of each country.

function profit_function(prim::Primitives, res::Results, X_prime::Int64, Q_value::Float64, ϵ_value::Float64)
    @unpack θ, α_N, α_K, l, r = prim
    @unpack C_star = res

    w = α_N *((θ-1)/θ) * l^(α_N *(θ-1)/θ-1)
    X = (1 + X_prime * (exp(Q_value)^θ) * C_star)^(1/θ) * exp(ϵ_value)

    a = α_N * (θ - 1)/θ
    b = α_K * (θ - 1)/θ
    Π = (1 - a - b) * X^(1/(1-a-b)) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b));

    Π
end

## calculate the profit of each country.

function Bellman(prim::Primitives, res::Results)
    @unpack r, T, N, Π_Q, n_Q, Q_grid, R = prim
    @unpack val_func, pol_func, Π_ϵ, n_ϵ, ϵ_grid = res



    val_func_new = zeros(2, n_ϵ, n_Q)
    pol_func_new = zeros(Int64, 2, n_ϵ, n_Q)
    median_rev = plant_revenue(prim, res, 0, 0.0, 0.0)[1] # non-exporting, ln(ϵ) = 0, ln(Q) = 0
    for ϵ_index = 1:n_ϵ
        for Q_index = 1:n_Q
            Π_0 = profit_function(prim, res, 0, Q_grid[Q_index], ϵ_grid[ϵ_index])
            Π_1 = profit_function(prim, res, 1, Q_grid[Q_index], ϵ_grid[ϵ_index])

            for X_index = 1:2
                val_func_new_ex = Π_1 - fixed_cost(prim, res, X_index-1, 1) * median_rev
                val_func_new_nx = Π_0
                for ϵp_index = 1:n_ϵ
                    for Qp_index = 1:n_Q
                        val_func_new_ex += 1 / (1 + R) * Π_ϵ[ϵ_index, ϵp_index] * Π_Q[Q_index, Qp_index] * val_func[2, ϵp_index, Qp_index]
                        val_func_new_nx += 1 / (1 + R) * Π_ϵ[ϵ_index, ϵp_index] * Π_Q[Q_index, Qp_index] * val_func[1, ϵp_index, Qp_index]
                    end
                end

                if val_func_new_ex > val_func_new_nx
                    pol_func_new[X_index, ϵ_index, Q_index] = 1
                    val_func_new[X_index, ϵ_index, Q_index] = val_func_new_ex
                else
                    pol_func_new[X_index, ϵ_index, Q_index] = 0
                    val_func_new[X_index, ϵ_index, Q_index] = val_func_new_nx
                end
            end
        end
    end

    val_func_new, pol_func_new
end

## define find_wage function to minimize the excess trade flow

function solve_model(prim::Primitives, res::Results)
    tol = 1e-4
    err = 1.0
    n = 1 # count

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

##

function get_index(val::Float64, grid::Array{Float64, 1})

    interp =  interpolate(grid, BSpline(Linear()))
    find_index(k) = abs(interp(k) - val)
    index = optimize(find_index, 1.0, length(grid)).minimizer
    index = round(Int, index)

    index
end

##

function simulate_Q(prim::Primitives, res::Results)
    @unpack T, ρ_Q, σ_Q, drop = prim

    Random.seed!(12341234)
    dist = Normal(0, σ_Q)
    Q = zeros(T)

    η = randn() * σ_Q/sqrt(1-ρ_Q)
    Q[1] = η

    for T_index = 2:T
        Q[T_index] = ρ_Q*Q[T_index-1] + randn() * σ_Q
    end

    Q
end

##

function simulate_epsilon(prim::Primitives, res::Results)
    @unpack T, N, drop = prim
    @unpack ρ_ϵ, σ_ϵ = res

    Random.seed!(12341234)
    dist = Normal(0, σ_ϵ)
    # η = randn()
    ϵ = zeros(T, N)

    # ϵ[1, :] = η[1, :]

    for N_index = 1:N
        η = randn() * σ_ϵ/sqrt(1-ρ_ϵ)
        ϵ[1, N_index] = η
        for T_index = 2:T
            ϵ[T_index, N_index] = ρ_ϵ*ϵ[T_index-1, N_index] + randn() * σ_ϵ
        end
    end

    ϵ
end

##

function simulate_pol_func(prim::Primitives, res::Results)
    @unpack N, T, Q_grid = prim
    @unpack pol_func, ϵ_grid = res
    @unpack ϵ_sim, Q_sim = sim

    pol_func_sim = zeros(Int64, T, N);
    domestic = zeros(T, N);
    exports = zeros(T, N);

    for N_index = 1:N
        for T_index = 1:T
            ϵ_index = get_index(ϵ_sim[T_index, N_index], ϵ_grid)
            ϵ_value = ϵ_grid[ϵ_index]
            Q_index = get_index(Q_sim[T_index], Q_grid)
            Q_value = Q_grid[Q_index]

            if T_index == 1
                pol_func_sim[T_index, N_index] = 0
            else
                pol_func_sim[T_index, N_index] = pol_func[pol_func_sim[T_index-1, N_index]+1, ϵ_index, Q_index]
            end
            domestic[T_index, N_index] = plant_revenue(prim, res, pol_func_sim[T_index, N_index], Q_value, ϵ_value)[2]
            exports[T_index, N_index] = plant_revenue(prim, res, pol_func_sim[T_index, N_index], Q_value, ϵ_value)[3]
        end
    end

    pol_func_sim, domestic, exports
end

##

function simulate(prim::Primitives, res::Results)
    @unpack T, N, drop = prim

    sim.ϵ_sim = simulate_epsilon(prim, res)
    sim.Q_sim = simulate_Q(prim, res)

    pol_func_sim, domestic_sim, exports_sim = simulate_pol_func(prim, res)

    starterrate = zeros(T - drop);
    stopperrate = zeros(T - drop);

    X = pol_func_sim[drop, :]
    pol_func_sim = pol_func_sim[drop+1:T, :]
    domestic_sim = domestic_sim[drop+1:T, :]
    exports_sim = exports_sim[drop+1:T, :]

    sim.pol_func_sim = pol_func_sim
    sim.domestic_sim = domestic_sim
    sim.exports_sim = exports_sim

    for T_index = 1:(T - drop)
        X_prime = pol_func_sim[T_index, :]
        starter = 0
        stopper = 0
        total_nonexperter = 1914 - sum(pol_func_sim[T_index, :])

        for N_index = 1:prim.N
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
        starterrate[T_index] = starter/total_nonexperter
        stopperrate[T_index] = stopper/(N - total_nonexperter)
        X = X_prime
    end

    starterrate = mean(starterrate)
    stopperrate = mean(stopperrate)

    Y = floor(Int64, (T-drop)/4)

    domestic_sim_annual = zeros(Y, N)
    exports_sim_annual = zeros(Y, N)
    exportsales = zeros(Y, N)

    for year = 1:Y
        for quarter = 1:4
            domestic_sim_annual[year, :] += domestic_sim[4*(year - 1) + quarter, :]
            exports_sim_annual[year, :] += exports_sim[4*(year - 1) + quarter, :]
        end
    end

    for N_index = 1:N
        for T_index = 1:Y
        exportsales[T_index, N_index] = exports_sim_annual[T_index, N_index]/(domestic_sim_annual[T_index, N_index] + exports_sim_annual[T_index, N_index])
        end
    end

    exporter = (exportsales .> 0)
    exportsales = sum(exportsales)/sum(exporter)

    cv = sqrt(var(domestic_sim_annual))/mean(domestic_sim_annual)

    y = zeros(N * Y)
    y_lag = zeros(N * Y)
    firmno = zeros(N * Y)
    year = zeros(N * Y)

    for N_index = 1:N
        for Y_index = 1:Y
            y[(N_index-1) * 100 + Y_index] = domestic_sim_annual[Y_index, N_index]
            firmno[(N_index-1) * 100 + Y_index] = N_index
            year[(N_index-1) * 100 + Y_index] = Y_index
        end
        for Y_index = 2:Y
            y_lag[(N_index-1) * 100 + Y_index] = domestic_sim_annual[Y_index-1, N_index]
        end
        y_lag[(N_index-1) * 100 + 1] = NaN
    end

    df = DataFrame(y_lag=y_lag, y=y, firmno=firmno, year=year)
    df = filter(row -> ! isnan(row.y_lag), df)

    γ = reg(df, @formula(y ~ y_lag + fe(firmno) + fe(year))).coef[1]

    [starterrate, stopperrate, exportsales, cv, γ]
end

##

function m_hat(f₀, f₁, C_star, σ_ϵ, ρ_ϵ)
     res.f₀ = f₀
     res.f₁ = f₁
     res.C_star = C_star
     res.σ_ϵ = σ_ϵ
     res.ρ_ϵ = ρ_ϵ

     solve_model(prim, res)
     simulate(prim, res)
end

##

function MoM(prim::Primitives, res::Results, sim::Simulations)
    @unpack f₀, f₁, C_star, σ_ϵ, ρ_ϵ = res
    @unpack moments_sim = sim

    W = Matrix{Float64}(I, 5, 5)
    J(f₀, f₁, C_star, σ_ϵ, ρ_ϵ) = (prim.moments_data - m_hat(f₀, f₁, C_star, σ_ϵ, ρ_ϵ))' * W * (prim.moments_data - m_hat(f₀, f₁, C_star, σ_ϵ, ρ_ϵ))
    lower = [0.1, 0.0, 0.0, 0.0, 0.3]
    upper = [1.0, 0.5, 0.6, 0.3, 1.0]
    # b_initial = (lower + upper)./2
    b_hat = optimize(b -> J(b[1], b[2], b[3], b[4], b[5]), lower, upper, [f₀, f₁, C_star, σ_ϵ, ρ_ϵ]).minimizer

    b_hat
end

##

MoM(prim, res, sim)
