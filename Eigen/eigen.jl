# hydrogen eigen solver

using SparseArrays, ArnoldiMethod

coulumb(x, y, z) = -1 / (4π * hypot(x, y, z))

function laplacian()
    I = [1:d;
         2:d; [1];
         res+1:d; 1:res;
         res^2+1:d; 1:res^2;
         1:d;
         1:d;
         1:d]

    J = [1:d;
         1:d;
         1:d;
         1:d;
         2:d; [1];
         res+1:d; 1:res;
         res^2+1:d; 1:res^2]

    V = [[-6 for i in 1:d];
         [1 for i in 1:6d]]

    if !periodic
        V[2d] = 0
        V[5d] = 0
    end

    sparse(I, J, V)
end

function potential()
    I = [1:d...]
    J = [1:d...]
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
    evec = open("data/evec$(j)$(p)_res$(res).dat", "w")
    for (x, y, z, A) in qs
        println(evec, "c3 $x $y $z $(A*scaler)")
    end
    println(evec, "T -0.9 0.9")
    println(evec, "E[$j] = $(E[j])")
    println(evec, "F")
end


const res = 100
const d = res^3
const O = (res + 1)/2

const qmax = 1e3
const dq = 2 * qmax / res
const λ = 1 / (2dq^2)

const thresh = 5e-5
const scaler = 1e3

const periodic = true
const p = periodic ? "p" : ""

const ψn = parse(Int64, ARGS[1])

function main()
    println("\ncalculating hamiltonian\n")

    Δ = laplacian()
    V = potential()

    H = -λ * Δ + V

    println("solving Hψ = Eψ\n")

    decomp, history = partialschur(H, nev=ψn, tol=1e-6, which=SR())

    E, ψ = partialeigen(decomp)

    println("solved!\n")

    Emin = minimum(E)

    println("Emin = $Emin\n")

    println("writing data\n")

    # enrdat = open("$(pwd())/eigen_data/energy$(periodic ? "p" : "").dat")

    # for i in eachindex(E) println(enrdat, i, " ", E[i]) end

    for j in 1:size(ψ, 2) anim(ψ[:,j], E, j) end

    println("finished in:")
end

@time main()

println()
