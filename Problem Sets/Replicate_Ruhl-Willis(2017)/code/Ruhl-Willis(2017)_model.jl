## keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    r::Float64 = 0.109  # annual interest rate
    R::Float64 = (1 + r)^(1 / 4) - 1  # quarterly interest rate
    ρ_Q::Float64 = 0.826  # real effective exchange rate process
    σ_Q::Float64 = 0.036  # real effective exchange rate process
    n_Q::Int64 = 11
    Π_Q::Array{Float64,2} = tauchen(n_Q, ρ_Q, σ_Q, 0.0)[1]
    Q_grid::Array{Float64} = tauchen(n_Q, ρ_Q, σ_Q, 0.0)[2]
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
    ϕ::Array{Float64}
    n_ϵ::Int64
    Π_ϵ::Array{Float64,2}
    ϵ_grid::Array{Float64}
    val_func::Array{Float64,3}
    pol_func::Array{Int64,3}
end


## structure that holds model results

mutable struct Simulations
    Q_sim::Array{Float64}
    ϵ_sim::Array{Float64, 2}
    moments_sim::Array{Float64}
end


## function for initializing model primitives and results

function Initialize()
    prim = Primitives()
    ϕ = [0.961, 0.04417998355827739, 0.14499085812207024, 0.14261294715764294, 0.862488261641639]
    ϕ = [0.961, 0.047, 0.146, 0.117, 0.873]
    # ϕ = [0.961, 0.047, 0.1440402171160336, 0.12833318683668238, 0.8600259464524087]    
    n_ϵ = 455
    Π_ϵ, ϵ_grid = tauchen(n_ϵ, ϕ[5], ϕ[4], 0.0)
    val_func = zeros(2, n_ϵ, prim.n_Q)
    pol_func = zeros(Int64, 2, n_ϵ, prim.n_Q)
    Q_sim = simulate_Q(prim)
    ϵ_sim = zeros(prim.T, prim.N)
    moments_sim = zeros(length(prim.moments_data))
    res = Results(ϕ, n_ϵ, Π_ϵ, ϵ_grid, val_func, pol_func)
    sim = Simulations(Q_sim, ϵ_sim, moments_sim)
    prim, res, sim
end


##

function tauchen(N::Int64, ρ::Float64, σ::Float64, μ::Float64)
    # Get discretized space using Tauchen's method (https://julia.quantecon.org/tools_and_techniques/finite_markov.html)

    dist = Normal(0, σ)
    σ_y = sqrt(σ^2 / (1 - ρ^2))

    y_bar = 3 * σ_y
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
    y = y .+ μ / (1 - ρ)
    Π, y
end


## calculate the cutoff productivity of each country pair

function fixed_cost(prim::Primitives, res::Results, X::Int64, X_prime::Int64)
    @unpack ϕ = res

    if X_prime == 0
        f = 0.0
    elseif X_prime == 1
        if X == 0
            f = ϕ[1]
        elseif X == 1
            f = ϕ[2]
        end
    end

    f
end


##

function plant_revenue(prim::Primitives, res::Results, X_prime::Int64, Q_value::Float64, ϵ_value::Float64)
    @unpack θ, α_N, α_K, l, r = prim
    @unpack ϕ = res

    C_star = ϕ[3]
    w = α_N * ((θ - 1) / θ) * l^(α_N * (θ - 1) / θ - 1)
    X = (1 + X_prime * (exp(Q_value)^θ) * C_star)^(1 / θ) * exp(ϵ_value)

    a = α_N * (θ - 1) / θ
    b = α_K * (θ - 1) / θ
    Y = exp(ϵ_value) * X^((a + b) / (1 - a - b)) * (a / w)^(a / (1 - a - b)) * (b / r)^(b / (1 - a - b))
    domestic = (1 / (X_prime * (exp(Q_value)^θ) * C_star + 1))^((θ - 1) / θ) * Y
    exports = X_prime * exp(Q_value) * (C_star^(1 / θ)) * (1 / (1 + exp(Q_value)^(-θ) / C_star))^((θ - 1) / θ) * Y

    Y, domestic, exports
end


## calculate the price index of each country.

function profit_function(prim::Primitives, res::Results, X_prime::Int64, Q_value::Float64, ϵ_value::Float64)
    @unpack θ, α_N, α_K, l, r = prim
    @unpack ϕ = res

    C_star = ϕ[3]
    w = α_N * ((θ - 1) / θ) * l^(α_N * (θ - 1) / θ - 1)
    X = (1 + X_prime * (exp(Q_value)^θ) * C_star)^(1 / θ) * exp(ϵ_value)

    a = α_N * (θ - 1) / θ
    b = α_K * (θ - 1) / θ
    Π = (1 - a - b) * X^(1 / (1 - a - b)) * (a / w)^(a / (1 - a - b)) * (b / r)^(b / (1 - a - b))

    Π
end


## calculate the profit of each country.

function Bellman(prim::Primitives, res::Results, Π_0, Π_1)
    @unpack r, T, N, Π_Q, n_Q, Q_grid, R = prim
    @unpack val_func, pol_func, Π_ϵ, n_ϵ, ϵ_grid = res

    val_func_new = zeros(2, n_ϵ, n_Q)
    pol_func_new = zeros(Int64, 2, n_ϵ, n_Q)
    # median_rev = plant_revenue(prim, res, 0, 0.0, 0.0)[1] # non-exporting, ln(ϵ) = 0, ln(Q) = 0

    @threads for ϵ_index = 1:n_ϵ
        for Q_index = 1:n_Q
            for X_index = 1:2
                val_func_new_ex = Π_1[Q_index, ϵ_index, X_index]
                val_func_new_nx = Π_0[Q_index, ϵ_index]

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
    @unpack n_Q, Q_grid = prim
    @unpack val_func, pol_func, Π_ϵ, n_ϵ, ϵ_grid = res

    tol = 1e-2
    err = 1.0
    n = 1 # count

    Π_0 = zeros(n_Q, n_ϵ)
    Π_1 = zeros(n_Q, n_ϵ, 2)
    median_rev = plant_revenue(prim, res, 0, 0.0, 0.0)[1] # non-exporting, ln(ϵ) = 0, ln(Q) = 0

    @threads for ϵ_index = 1:n_ϵ
        for Q_index = 1:n_Q
            Π_0[Q_index, ϵ_index] = profit_function(prim, res, 0, Q_grid[Q_index], ϵ_grid[ϵ_index])
            for X_index = 1:2
                Π_1[Q_index, ϵ_index, X_index] = profit_function(prim, res, 1, Q_grid[Q_index], ϵ_grid[ϵ_index]) - fixed_cost(prim, res, X_index - 1, 1) * median_rev

            end
        end
    end

    while true
        val_func_new, pol_func_new = Bellman(prim, res, Π_0, Π_1)
        err = maximum(abs.(res.val_func - val_func_new))
        res.val_func = val_func_new
        res.pol_func = pol_func_new
        # println("***** ", n, "th iteration *****")
        # @printf("Absolute difference: %0.5f \n", float(err))
        # println("***************************")
        n += 1
        if err < tol
            break
        end
    end
    println("Value function converged after ", n, " iterations.")
    # println(" ")

end

##

function get_index(val::Float64, grid::Array{Float64,1})

    interp = interpolate(grid, BSpline(Linear()))
    find_index(k) = abs(interp(k) - val)
    index = optimize(find_index, 1.0, length(grid)).minimizer
    index = round(Int, index)

    index
end

##

function simulate_Q(prim::Primitives)
    @unpack T, ρ_Q, σ_Q, drop = prim

    Random.seed!(12341234)
    # dist = Normal(0, σ_Q)
    Q = zeros(T)

    η = randn() * σ_Q / sqrt(1 - ρ_Q)
    Q[1] = η

    for T_index = 2:T
        Q[T_index] = ρ_Q * Q[T_index-1] + randn() * σ_Q
    end

    Q
end


function simulate_epsilon(prim::Primitives, res::Results)
    @unpack T, N, drop = prim
    @unpack ϕ = res
    
    σ_ϵ, ρ_ϵ = ϕ[4], ϕ[5]

    Random.seed!(12)
 
    ϵ = zeros(T, N)
 
    η = randn(N) * σ_ϵ / sqrt(1 - ρ_ϵ)
    ϵ[1, :] = η
    for T_index = 2:T
        ϵ[T_index, :] = ρ_ϵ * ϵ[T_index-1, :] .+ randn(N) * σ_ϵ
    end


    ϵ
end


##

function simulate_pol_func(prim::Primitives, res::Results)
    @unpack N, T, Q_grid = prim
    @unpack pol_func, ϵ_grid = res
    @unpack ϵ_sim, Q_sim = sim

    pol_func_sim = zeros(Int64, T, N)
    domestic = zeros(T, N)
    exports = zeros(T, N)

    @threads for N_index = 1:N
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
    Y = floor(Int64, (T - drop) / 4)

    sim.ϵ_sim = simulate_epsilon(prim, res)
    # sim.ϵ_sim = zeros(T, N)
    # for N_index = 1:N
    #     sim.ϵ_sim[:, N_index] = AR1sim(T, res.ϵ_grid, res.Π_ϵ)
    # end

    pol_func_sim_temp, domestic_sim_temp, exports_sim_temp = simulate_pol_func(prim, res)

    # starterrate = zeros(T - drop)
    # stopperrate = zeros(T - drop)

    # X = pol_func_sim[drop, :]
    # for T_index = 1:(T-drop)
    #     X = pol_func_sim[T_index, :]
    #     X_prime = pol_func_sim[T_index+1, :]
    #     starter = 0
    #     stopper = 0
    #     total_nonexperter = 1914 - sum(X)

    #     for N_index = 1:prim.N
    #         if X[N_index] == 0
    #             if X_prime[N_index] == 1
    #                 starter = starter + 1
    #             end
    #         else
    #             if X_prime[N_index] == 0
    #                 stopper = stopper + 1
    #             end
    #         end
    #     end
    #     starterrate[T_index] = starter / total_nonexperter
    #     stopperrate[T_index] = stopper / (N - total_nonexperter)
    # end

    pol_func_sim = pol_func_sim_temp[drop+1:T, :]
    domestic_sim = domestic_sim_temp[drop+1:T, :]
    exports_sim = exports_sim_temp[drop+1:T, :]

    pol_func_sim_annual = zeros(Y, N)
    domestic_sim_annual = zeros(Y, N)
    exports_sim_annual = zeros(Y, N)

    @threads for year = 1:Y
        for quarter = 1:4
            pol_func_sim_annual[year, :] += pol_func_sim[4*(year-1)+quarter, :]
            domestic_sim_annual[year, :] += domestic_sim[4*(year-1)+quarter, :]
            exports_sim_annual[year, :] += exports_sim[4*(year-1)+quarter, :]
        end
    end
    pol_func_sim_annual = pol_func_sim_annual .> 0

    starterrate = zeros(Y - 1)
    stopperrate = zeros(Y - 1)

    @threads for Y_index = 1:Y-1
        X = pol_func_sim_annual[Y_index, :]
        X_prime = pol_func_sim_annual[Y_index+1, :]
        starter = 0
        stopper = 0
        total_nonexporter = N - sum(X)

        if total_nonexporter == 0 || total_nonexporter == N
            starterrate[Y_index] = 0.0
            stopperrate[Y_index] = 0.0
        else
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
            starterrate[Y_index] = starter / total_nonexporter
            stopperrate[Y_index] = stopper / (N - total_nonexporter)
        end
    end
    m1 = mean(starterrate)
    m2 = mean(stopperrate)

    exportsales = exports_sim_annual ./ (domestic_sim_annual + exports_sim_annual)
    exporter = (exportsales .> 0)
    m3 = sum(exportsales) / sum(exporter)
    if sum(exporter) == 0
        m3 = 0.0
    end
    # exportsales = sum(exportsales) / sum(exporter)
a = domestic_sim + exports_sim
    cv = zeros(400)
    cv = var(log.(a[:, :])) / mean(log.(a[:, :]))

    @threads for Y_index = 1:400
        cv[Y_index] = var(log.(a[Y_index, :])) / mean(log.(a[Y_index, :]))
    end
    log(1.4)
log(1.6)
    std(domestic_sim)/mean(domestic_sim)
    m4 = mean(cv)
    # cv = zeros(400)
    # for T_index = 1:400
    #     cv[T_index] = std(log.(domestic_sim[T_index, :])) / mean(log.(domestic_sim[T_index, :]))
    # end

    y = zeros(N * Y)
    y_lag = zeros(N * Y)
    firmno = zeros(N * Y)
    year = zeros(N * Y)

    @threads for N_index = 1:N
        for Y_index = 1:Y
            y[(N_index-1)*100+Y_index] = log.(domestic_sim_annual[Y_index, N_index])
            firmno[(N_index-1)*100+Y_index] = N_index
            year[(N_index-1)*100+Y_index] = Y_index
        end
        for Y_index = 2:Y
            y_lag[(N_index-1)*100+Y_index] = log.(domestic_sim_annual[Y_index-1, N_index])
        end
        y_lag[(N_index-1)*100+1] = NaN
    end


    df = DataFrame(y_lag=y_lag, y=y, firmno=firmno, year=year)
    df = filter(row -> !isnan(row.y_lag), df)

    γ = reg(df, @formula(y ~ y_lag + fe(firmno) + fe(year))).coef[1]

    [m1, m2, m3, m4, γ]
end

##
function objective_function(prim, param)
    @unpack moments_data = prim

    println(param)

    if param[1] < 0 || param[1] > 1 || param[2] < 0 || param[2] > 1 ||param[3] < 0 || param[4] < 0 || param[5] < 0 || param[5] > 1 
        ObjFn = 10e100

    else
        res.ϕ = param

        # W = Matrix{Float64}(I, 5, 5)

        solve_model(prim, res)
        M = simulate(prim, res)

        ObjFn = (moments_data - M)' * (moments_data - M)
        
        sim.moments_sim = M
    end

    return ObjFn
end

##

function Minimize_Moment(prim::Primitives, res::Results, initial)

    opt = optimize(phi -> objective_function(prim, phi), initial, BFGS(), Optim.Options(show_trace=true, iterations=10))

    opt
end
