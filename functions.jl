using JuMP
using MosekTools
import Dualization
using LinearAlgebra
import JLD2
using Ket

# Genera un vector de d + 1 bases MUB para 2 ≤ d ≤ 13. Para d = 6, 10, 12, las bases son solo aproximadamente no sesgadas
function numerical_mub(d::Int)
    mub_dict = JLD2.load("mubs.jld2")   
    return mub_dict["mubs"][d]
end

# Genera un vector de d + 1 proyectores asociados a que Alice y Bob obtengan el mismo resultado, construidos con las bases MUB
function E(d)
    E_matrices = [zeros(ComplexF64, d^2, d^2) for _ = 1:d+1]
    mub = numerical_mub(d)
    for k = 1:d+1
        for i = 1:d
            P = kron(ketbra(mub[k][:, i]), transpose(ketbra(mub[k][:, i])))
            E_matrices[k] += P
        end
    end
    cleanup!.(E_matrices)
    return E_matrices
end

# Estado isótropo 
function isotropo(d, v::Real)
    phi = zeros(d^2)
    for i = 1:d
        phi_i = 1 / sqrt(d) * kron(ket(i, d), ket(i, d))
        phi += phi_i
    end
    return v * ketbra(phi) + (1 - v) / d^2 * Matrix(I, d^2, d^2)
end

# Genera un vector de d + 1 probabilidades de que Alice y Bob obtengan resultados iguales
function f(d, v)
    prob = zeros(d + 1)
    for k = 1:d+1
        prob[k] = dot(isotropo(d, v), E(d)[k])
    end
    return prob
end

function entropia_binaria(p::Real)
    function log(x)
        if x == 0
            sol = 0.0
        else
            sol = log2(x)
        end
    end
    return -p * log(p) - (1 - p) * log(1 - p)
end

function gauss_radau(m::Int)
    diagonal = fill(0.5, m)
    diagonal[m] = (3 * m - 1) / (4 * m - 2)
    superdiagonal = [n / (2 * sqrt(4 * n^2 - 1)) for n = 1:m-1]
    J = SymTridiagonal(diagonal, superdiagonal)
    t, vecs = eigen(J)
    w = vecs[1, :] .^ 2
    return w, t
end

# Entropía condicional S(Alice/Eva) (límite inferior de la entropía relativa de von Neumann)
function Sae(m, d, v)
    model = Model()
    set_optimizer(model, Dualization.dual_optimizer(Mosek.Optimizer))
    #   set_optimizer(model, Mosek.Optimizer)
    
    w, t = gauss_radau(m)
    c_m = sum(w ./ (t .* log(2)))

    I_B = Matrix(I, d, d)

    @variable(model, sigma[1:d^2, 1:d^2], Hermitian)
    ζ = [@variable(model, [1:d^2, 1:d^2] in ComplexPlane()) for i = 1:m, a = 1:d]
    η = [@variable(model, [1:d^2, 1:d^2], Hermitian) for i = 1:m, a = 1:d]
    θ = [@variable(model, [1:d^2, 1:d^2], Hermitian) for i = 1:m, a = 1:d]

    @constraint(model, tr(sigma) == 1)
    for k = 1:d+1
        @constraint(model, dot(sigma, E(d)[k]) == f(d, v)[k])
    end
    for i = 1:m
        for a = 1:d
            @constraint(model, Hermitian([sigma ζ[i, a]; ζ[i, a]' η[i, a]]) in HermitianPSDCone())
            @constraint(model, Hermitian([sigma ζ[i, a]'; ζ[i, a] θ[i, a]]) in HermitianPSDCone())
        end
    end

    obj = c_m + sum(
        w[i] / (t[i] * log(2)) * (
            dot(
                kron(ketbra(ket(a, d)), I_B),
                Hermitian(
                    ζ[i, a] + ζ[i, a]'
                ) + 
                (1 - t[i]) * η[i, a]
            ) +
            t[i] * tr(θ[i, a])
        ) 
        for i = 1:m, a = 1:d
    )

    @objective(model, Min, real(obj))

    optimize!(model)
    return objective_value(model)
end

# Entropía condicional de Shannon H(A/B)
function Hab(d, v)
    if v == 1 || v == 0
        return 0.0
    else
        return entropia_binaria(v + (1 - v) / d) + (1 - v - (1 - v) / d) * log2(d - 1)
    end
end

# Tasa de rondas válidas de la implementación con grados de libertad temporales
function R_time(d, t_b, T, gamma)
    if T == Inf
        return 0.0
    else
        return (d * t_b * T^2 + gamma) * exp(- d * t_b * (2 * T + gamma))
    end
end

# Tasa de clave de la implementación con grados de libertad temporales
function K_time(m, d, noise_signal)
    t_b = 1.31 * 10^-9
    lambda = 2 * 10^5
    PLA = 0.01
    PLB = 0.984
    PCA = 0.6
    PCB = 0.6
    gamma = lambda * PCA * PCB * (1 - PLA) * (1 - PLB)   
    
    T = noise_signal * gamma / (1 - noise_signal)
    v = 1 / (1 + d * t_b * T^2 * gamma^-1) 
    
    return R_time(d, t_b, T, gamma) * (Sae(m, d, v) - Hab(d, v))
end

# Tasa de clave calculada con el método de min-entropía
# Acotamos inferiormente la entropía condicional de von Neumann con la min-entropía condicional
function K_time_min_entropy(d, s, noise_signal)
    t_b = 1.31 * 10^-9
    lambda = 2 * 10^5
    PLA = 0.01
    PLB = 0.984
    PCA = 0.6
    PCB = 0.6
    gamma = lambda * PCA * PCB * (1 - PLA) * (1 - PLB)   

    T = noise_signal * gamma / (1 - noise_signal)
    v = 1 / (1 + d * t_b * T^2 * gamma^-1) 
    
    Hae_min = -log2((sqrt(v * d + 1 - v) + (d - 1) * sqrt(1 -v))^2 / d^2)
    return R_time(d, t_b, T, gamma) * (Hae_min - Hab(d, v))
end


# Tasa de rondas válidas de la implementación con grados de libertad espaciales
function R_space(d, Δt, μ, ξ, gamma)
    return (d / Δt) * exp(- Δt * ((d - 1) * (2 * μ + 2 * ξ / d) + gamma)) *
    (d * (1 - exp(- Δt * (μ + ξ / d)))^2 + exp(Δt * gamma / d) - 1)
end

# Tasa de clave de la implementación con grados de libertad espaciales
function K_space(m, d, noise_signal)
    Δt = 10^-7
    PP = 1
    lambda = 2 * 10^5
    PLA = 0.01
    PLB = 0.984
    PCA = 0.6
    PCB = 0.6
    μ = 500
    gamma = PP * lambda * PCA * PCB * (1 - PLA) * (1 - PLB)   

    ξ = (μ * (noise_signal - 1) + noise_signal * gamma) / (1 - noise_signal)
    v = (exp(Δt * gamma / d) - 1) / (exp(Δt * gamma / d) - 1 + d * (1 - exp(- Δt * (μ + ξ / d)))^2) 

    return R_space(d, Δt, μ, ξ, gamma) * (Sae(m, d, v) - Hab(d, v))
end

function K_space_min_entropy(d, s, noise_signal)
    Δt = 10^-7
    PP = 1
    lambda = 2 * 10^5
    PLA = 0.01
    PLB = 0.984
    PCA = 0.6
    PCB = 0.6
    μ = 500
    gamma = PP * lambda * PCA * PCB * (1 - PLA) * (1 - PLB)   
    
    ξ = (μ * (noise_signal - 1) + noise_signal * gamma) / (1 - noise_signal)
    v = (exp(Δt * gamma / d) - 1) / (exp(Δt * gamma / d) - 1 + d * (1 - exp(- Δt * (μ + ξ / d)))^2) 

    K_iso_tot = (log2(s) - 2 * log2(sqrt(v * d + 1 - v) + 
    (s - 1)* sqrt(1 - v))) + (v * d + 1 - v) / (v * d + (1 - v) * s) * (log2(v * d + 1 - v)) + 
    (s - 1) * (1 - v) / (v * d + (1 - v) * s) * log2(1 - v)

    K_space_min_entropy = R_space(d, Δt, μ, ξ, gamma) * K_iso_tot 

    return K_space_min_entropy
end
