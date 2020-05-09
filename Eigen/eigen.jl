# hydrogen eigen solver

using ArnoldiMethod, LinearAlgebra, SparseArrays

coulumb(x, y, z) = -1 / (4π * hypot(x, y, z))

function laplacian()
    I = [1:dim;
         2:dim; # [1];
         res+1:dim; 1:res;
         res^2+1:dim; 1:res^2;
         1:dim-1;
         1:dim;
         1:dim]

    J = [1:dim;
         1:dim-1;
         1:dim;
         1:dim;
         2:dim; # [1];
         res+1:dim; 1:res;
         res^2+1:dim; 1:res^2]

    V = [[-6 for i in 1:dim];
         [1 for i in 1:6dim-2]]

    sparse(I, J, V)
end

function potential()
    I = [1:dim;]
    J = [1:dim;]
    V = Vector{Float64}(undef, dim)
    for k = 1:res
        for j = 1:res
            for i = 1:res
                x = i - origin
                y = j - origin
                z = k - origin
                V[i + (j - 1)*res + (k - 1)*res^2] = coulumb(x, y, z)
            end
        end
    end
    sparse(I, J, V)
end

function lattice(ψ::Vector)
    qs = []
    for k = 1:res
        for j = 1:res
            for i = 1:res
                A = abs2(ψ[i + (j - 1)*res + (k - 1)*res^2])
                x = i - origin
                y = j - origin
                z = k - origin
                if A > B / res push!(qs, (x, y, z, A)) end
            end
        end
    end
    qs
end

function anim(ψ::Vector, E::Float64, j::Int)
    qs = lattice(ψ)
    wfn = open("data/res$(res)_wfn$j.txt", "w")
    for (x, y, z, A) in qs
        println(wfn, "c3 $x $y $z $(A*R)")
    end
    println(wfn, "T -0.9 0.9")
    println(wfn, "E[$j] = $E")
    println(wfn, "F")
end

anim(E::Float64, j::Int) = ψ -> anim(ψ, E, j)

const res = parse(Int, ARGS[1])
const dim = res^3
const origin = (res - 1) / 2

const N = parse(Int, ARGS[2])

# A = abs2(ψ[i,j,k]) < B / res => do not animate ψ[i,j,k]
const B = 1e-3

# scale A by R before writing
const R = 3e3

function main()
    println("\ncalculating H...\n")

    Δ = laplacian()
    V = potential()
    H = -0.5Δ + V

    println("solving: Hψ = Eψ\n")

    decomp = partialschur(H,
                          nev=N,
                          tol=1e-4,
                          which=SR())[1]

    println("got shur decomp...\n")

    E, Ψ = partialeigen(decomp)

    println("solved!\n")

    n = size(Ψ, 2)

    println("found $n / $N wavefunctions")

    Emin = minimum(E)

    println("Emin = $Emin\n")
    println("writing wavefunctions to data\n")

    for j = 1:n Ψ[:,j] |> normalize! |> anim(E[j], j) end

    println("finished in:\n")
end

@time main()

println()

