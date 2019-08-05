# hydrogen eigen solver

using LinearAlgebra, SparseArrays, Arpack

coulumb(x, y, z) = -e^2 / (4π * ϵ₀ * hypot(x, y, z))

function laplacian()
    I = [1:dim;
         2:dim; [1];
         res+1:dim; 1:res;
         res^2+1:dim; 1:res^2;
         1:dim;
         1:dim;
         1:dim]

    J = [1:dim;
         1:dim;
         1:dim;
         1:dim;
         2:dim; [1];
         res+1:dim; 1:res;
         res^2+1:dim; 1:res^2]

    V = [[-6 for i in 1:dim];
         [1 for i in 1:(2 * 3 * dim)]]

    if !periodic
        V[2dim] = 0
        V[5dim] = 0
    end

    sparse(I, J, V)
end

function potential()
    I = [1:dim...]
    J = [1:dim...]
    V::Vector{Float64} = []
    for k = 1:res
        for j = 1:res
            for i = 1:res
                x = (i - O)*dq
                y = (j - O)*dq
                z = (k - O)*dq
                push!(V, coulumb(x, y, z))
            end
        end
    end
    sparse(I, J, V)
end

function lattice(ψ::Vector{Float64})
    qs = []
    for k = 1:res
        for j = 1:res
            for i = 1:res
                x = (i - O)*dq
                y = (j - O)*dq
                z = (k - O)*dq
                A = abs(ψ[(i + (j - 1)*res + (k - 1)*res^2)])^2
                if A > thresh
                    push!(qs, (x, y, z, A))
                end
            end
        end
    end
    qs
end

function anim(ψ::Vector, E::Vector{Float64}, j::Int)
    qs = lattice(ψ)
    p = periodic ? "p" : ""
    evec = open("evec$(j)$(p).dat", "w")
    for i in 1:itr
        for (x, y, z, A) in qs
            println(evec, "c3 $x $y $z $(A*scaler)")
        end
        println(evec, "T -0.9 0.9")
        println(evec, "E[$j] = $(E[j])")
        println(evec, "F")
    end
    println(evec, "Q")
end

const c = 3e8
const ħ = 1 #6.5821195e-16 #1.05457e-34
const e = 1
const m = 1 #5.109989561e5 * c^2
const ϵ₀ = 1 #8.8541878e-12

const res = 30
const dim = res^3
const qmax = 250
const dq = 2 * qmax / res

const thresh = 2.0e-4
const scaler = 5.5e3

const itr = 1000

const O = (res + 1)/2
const λ = ħ^2 / (2m * dq^2)
const ψn = parse(Int64, ARGS[1])

const periodic = true

function main()
    println("calculating hamiltonian\n")

    ∇² = laplacian()
    V  = potential()

    H = -λ * ∇² + V

    println("solving Hψ = Eψ\n")

    E, ψ = eigs(H, nev=ψn, which=:SM)

    println("solved!\n")

    Emin = minimum(E)

    println("Emin = $Emin\n")

    enrdat = open("energy$(periodic ? "p" : "").dat", "w")

    for i in eachindex(E) println(enrdat, i, " ", E[i]) end

    for j in 1:size(ψ, 2) anim(ψ[:,j], E, j) end

    println("finished in:")
end

@time main()
