## keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    r::Float64 = 0.109  # interest rate
    R::Float64 = (1+r)^(1/4)-1
    ρ_Q::Float64 = 0.826  # real effective exchange rate process
    σ_Q::Float64 = 0.036  # real effective exchange rate process
    n_Q::Int64 = 11
    Π_Q::Array{Float64, 2} = tauchen(n_Q, ρ_Q, σ_Q)[1]
    Q_grid::Array{Float64} = tauchen(n_Q, ρ_Q, σ_Q)[2]
    α_N::Float64 = 0.45  # labor share
    α_K::Float64 = 0.55  # capital share
    l::Int64 = 63  # median plant employee
    θ::Float64 = 5.0  # elasticity of substitution
    T::Int64 = 420
    N::Int64 = 1914
    drop::Int64 = 20
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
    Q::Array{Float64}
    ϵ::Array{Float64, 2}
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
    Π_ϵ = tauchen(n_ϵ, ρ_ϵ, σ_ϵ)[1]
    ϵ_grid = tauchen(n_ϵ, ρ_ϵ, σ_ϵ)[2]
    val_func = zeros(2, n_ϵ, prim.n_Q)
    pol_func = zeros(Int64, 2, n_ϵ, prim.n_Q)
    Q = zeros(prim.T - prim.drop)
    ϵ = zeros(prim.T - prim.drop, prim.N)
    res = Results(f₀, f₁, C_star, σ_ϵ, ρ_ϵ, n_ϵ, Π_ϵ, ϵ_grid, val_func, pol_func, Q, ϵ)
    prim, res
end

##

std_norm_cdf(x::T) where {T <: Real} = 0.5 * erfc(-x/sqrt(2))
std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* erfc(-x./sqrt(2))

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Integer=3) where {T1 <: Real, T2 <: Real}
    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    # NOTE: I need to shift this vector after finding probabilities
    #       because when finding the probabilities I use a function
    #       std_norm_cdf that assumes its input argument is distributed
    #       N(0, 1). After adding the mean E[y] is no longer 0, so
    #       I would be passing elements with the wrong distribution.
    #
    #       It is ok to do after the fact because adding this constant to each
    #       term effectively shifts the entire distribution. Because the
    #       normal distribution is symmetric and we just care about relative
    #       distances between points, the probabilities will be the same.
    #
    #       I could have shifted it before, but then I would need to evaluate
    #       the cdf with a function that allows the distribution of input
    #       arguments to be [μ/(1 - ρ), 1] instead of [0, 1]

    yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable
    grid = collect(yy)
    # renormalize. In some test cases the rows sum to something that is 2e-15
    # away from 1.0, which caused problems in the MarkovChain constructor
    Π = Π./sum(Π, dims = 2)

    Π, grid
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


## calculate the price index of each country.

function profit_function(prim::Primitives, res::Results, X_prime::Int64, Q_value::Float64, ϵ_value::Float64)
    @unpack θ, α_N, α_K, l, r = prim
    @unpack C_star = res

    w = α_N *((θ-1)/θ) * l^(α_N *(θ-1)/θ-1)
    X = (1 + X_prime * (exp(Q_value)^θ) * C_star)^(1/θ) * exp(ϵ_value)

    # X = (1 .+ 1*(Q.^θ) * C_star).^(1/θ) .* ϵ[3,1]
    a = α_N * (θ - 1)/θ
    b = α_K * (θ - 1)/θ
    Π = (1 - a - b) * X^(1/(1-a-b)) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b));

    Π
end

## calculate the profit of each country.

function Bellman(prim::Primitives, res::Results)
    @unpack r, T, N, Π_Q, n_Q, Q_grid = prim
    @unpack val_func, pol_func, Π_ϵ, n_ϵ, ϵ_grid = res

    val_func_new = zeros(2, n_ϵ, n_Q)
    pol_func_new = zeros(Int64, 2, n_ϵ, n_Q)

    for ϵ_index = 1:n_ϵ
        for Q_index = 1:n_Q
            Π_0 = profit_function(prim, res, 0, Q_grid[Q_index], ϵ_grid[ϵ_index])
            Π_1 = profit_function(prim, res, 1, Q_grid[Q_index], ϵ_grid[ϵ_index])

            for X_index = 1:2
                value_export = Π_1 - Π_0 - fixed_cost(prim, res, X_index-1, 1)
                for ϵp_index = 1:n_ϵ, Qp_index = 1:n_Q
                    value_export += 1/(1+r) * Π_ϵ[ϵ_index, ϵp_index] * Π_Q[Q_index, Qp_index] * (val_func[2, ϵp_index, Qp_index] - val_func[1, ϵp_index, Qp_index])
                end
                if value_export>0
                    pol_func_new[X_index, ϵ_index, Q_index] = 1
                    val_func_new[X_index, ϵ_index, Q_index] = Π_1 - fixed_cost(prim, res, X_index-1, 1)
                else
                    val_func_new[X_index, ϵ_index, Q_index] = Π_0
                end

                for ϵp_index = 1:n_ϵ, Qp_index = 1:n_Q
                        val_func_new[X_index, ϵ_index, Q_index] += 1/(1+r) * Π_ϵ[ϵ_index, ϵp_index] * Π_Q[Q_index, Qp_index] * (val_func[pol_func_new[X_index, ϵ_index, Q_index]+1, ϵp_index, Qp_index])
                end
            end
        end

        # println(ϵ_index)
    end

    val_func_new, pol_func_new
end

## define find_wage function to minimize the excess trade flow

function solve_model(prim::Primitives, res::Results)
    # tol = 1e-4
    tol = 1
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

function Simulate_Q(prim::Primitives, res::Results)
    @unpack T, ρ_Q, σ_Q, drop = prim

    Random.seed!(12341234)
    dist = Normal(0, σ_Q)
    ϵ = rand(dist, T)
    Q = zeros(T)

    Q[1] = ϵ[1]

    for T_index = 2:T
        Q[T_index] = ρ_Q*Q[T_index-1] + ϵ[T_index]
    end

    plot_Q = plot(collect(1:T), Q, labels = "", title = "Movement of Q")

    Q = exp.(Q)

    res.Q = Q
end

##

function Simulate_epsilon(prim::Primitives, res::Results)
    @unpack T, N, drop = prim
    @unpack ρ_ϵ, σ_ϵ = res

    Random.seed!(12341234)
    dist = Normal(0, σ_ϵ)
    η = rand(dist, T, N)
    ϵ = zeros(T, N)

    ϵ[1, :] = η[1, :]

    for N_index = 1:N
        for T_index = 2:T
            ϵ[T_index, N_index] = ρ_ϵ*ϵ[T_index-1, N_index] + η[T_index, N_index]
        end
    end

    plot_ϵ = plot(collect(1:T), ϵ[:, 3], labels = "", title = "Movement of ϵ")

    ϵ = exp.(ϵ)

    res.ϵ = ϵ
end
