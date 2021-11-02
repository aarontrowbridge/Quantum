# hydrogen eigen solver

using Arpack, LinearAlgebra, SparseArrays

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

# define ⊗ operator as kron
⊗(A, B) = kron(A, B)

# create laplacian matrix from kron products of second derivative matricies
function kron_laplacian(res)
    ∂² = Tridiagonal(ones(res - 1), fill(-2.0, (res)), ones(res - 1))
    ∂
    Id = Matrix{Float64}(I, res, res)

    return (∂² ⊗ Id ⊗ Id) + (Id ⊗ ∂² ⊗ Id) + (Id ⊗ Id ⊗ ∂²)
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


# number of eigenvectors
# const N = parse(Int, ARGS[2])

# lattice resolution
# const res = parse(Int, ARGS[1])

const res = 3

const dim = res^3

const origin = (res - 1) / 2

# A = abs2(ψ[i,j,k]) < B / res => do not animate ψ[i,j,k]
const B = 1e-3

# scale A by R before writing
const R = 3e3



print(kron_laplacian(res))






# function main()
#     println("\ncalculating H...\n")

#     Δ = laplacian()
#     V = potential()
#     H = -0.5Δ + V

#     println("solving: Hψ = Eψ\n")

#     E, Ψ = eigs(H, howmany=N)
    
#     println("solved!\n")

#     Emin = minimum(E)

#     println("Emin = $Emin\n")
#     println("writing wavefunctions to data\n")

#     for j = 1:N Ψ[:,j] |> normalize! |> anim(E[j], j) end

#     println("finished in:\n")
# end

# @time main()

# println()

