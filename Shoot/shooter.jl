# shooting method
using Printf, LinearAlgebra

coulumb(r::Float64) = -1 / r + l*(l + 1) / r^2
quadratic(x::Float64) = (x - origin)^2

f(x, E) = 2(E - V(x))

function shoot(E::Float64)
    u = zeros(res)
    u[res] = 0.0
    u[res - 1] = 1.0
    q₀ = 0.0
    q₁ = 1 + dr2/12 * f(R - dr, E)
    for i = (res - 2):-1:1
        q₂ = 2q₁ - q₀ - dr2 * f((i+1)*dr, E) * u[i+1]
        q₀ = q₁
        q₁ = q₂
        u[i] = (q₂ / (1 + dr2/12 * f(i*dr, E)))
    end
    normalize!(u)
    10u
end

function bisect(lbnd::Float64, ubnd::Float64, n)
    a = lbnd
    b = ubnd
    itr = 1
    while itr <= max
        c = (a + b) / 2
        ψ₀ = shoot(c)
        printer(ψ₀, n, c)
        if abs(ψ₀[1])^2 < tol
            @printf "!hit tolerence!\n"
            @printf "!E = %f\n!ψ(0) = %e\n!\n" c ψ₀[1]
            return c
            break
        end
        itr += 1
        if itr == max
            @printf "!reached max iterations!\n"
            @printf "!E = %f\n!ψ(0) = %e\n!\n" c ψ₀[1]
            return c
            break
        end
        ψ₁ = shoot(a)
        if sign(ψ₀[1]) == sign(ψ₁[1])
            a = c
        else
            b = c
        end
    end
end

function printer(ψ::Vector, n, E)
    @printf "F\n"
    @printf "A 1\n"
    @printf "l 0.0 0.0 15.0 0.0\n"
    @printf "l 0.0 -15.0 0.0 15.0\n"
    @printf "T 0.1 0.7\n"
    @printf "n = %i\n" n
    @printf "T 0.1 0.6\n"
    @printf "E = %f\n" E
    for i = eachindex(ψ)
        @printf "c3 %f %f 0.0 0.03\n" i*dr ψ[i]
    end
end

const R = 15
const origin = R/2
const res = 1500
const dr = R / res
const dr2 = dr^2

const tol = 1.0e-8
const max = 200

const V = quadratic
const radial = false

const N = 50

const E₀ = -1.0
const dE = 0.1

function main()
    E₁ = E₀
    n = 1
    while n <= N
        E₂ = E₁ + dE
        ψ₁ = shoot(E₁)
        ψ₂ = shoot(E₂)
        printer(ψ₁, n, E₁)
        printer(ψ₂, n, E₂)
        if sign(ψ₁[1]) == sign(ψ₂[1])
            E₁ = E₂
        else
            @printf "!shooting! n = %i, E ∈ [%g, %g]\n" n E₁ E₂
            E = bisect(E₁, E₂, n)
            E₁ = E₂
            n += 1
        end
    end
    @printf "Q\n"
end

@time main()

