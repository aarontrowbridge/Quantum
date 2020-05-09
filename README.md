# Quantum Mechanical Simulations
A collection of scripts to find simulated position space eigen vectors, or wavefunctions, of the discretized time independent and dependent Sch&ouml;dinger Equation. 

The time Indpendent (Implicit) Schr&ouml;dinger Equation: *ISE*
```
ISE: H|psi(x)> = E|psi(x)> 
```

The time Dependent (Deterministic) Schr&ouml;dinger Equation: *DSE*
```
DSE: |psi(x, t)> = exp(im*H*t)|psi(x, t=0)>
```

## Shoot
The Shooting Method, or "wagging the tail", is an algorithm to numerically solve for the wavefunctions of the *ISE*, for arbitrary potentials in 1-D, including the radial equation for spherically symmetric potentials. Here's what it looks like in julia:

```julia
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
        u[i] = q₂ / (1 + dr2/12 * f(i*dr, E))
    end
    normalize!(u)
    return u
end
```

The idea is that the we treat the *ISE* as a boundary value problem, where *psi*(0) = 0 = *psi*(R); we then guess a low value for E, and actually integrate the *ISE* from *R*in to 0 , which is a "shot"; we then see if the *psi*(0) < or > 0; we then step, *E* <- *E* + *dE*, until we cross over 0, and then use bisection to narrow down a solution. When animated it looks like the wavefunction is wagging its tail, hence the popular name. 

A novel approach taken here is that instead of using Euler's method for numerical integration, this simulation implements Numerov's method, a 4th order implicit method suited perfectly for the *ISE*.

I need to figure out how to include `LaTeX` here.

## Eigen
This method "simply" sets up the matrix eigen value problem associated with solving the *ISE* on a 3-D lattice, where each dimension has a resolution of *res* intervals, resulting in a cube in space with *res*^3 lattice points.

This can be flattened into a vector of length *res*^3, with a corresponding Hamiltonian operator, *H* = 0.5*L* + *V*, where *L* is the 3-D Laplacian operator and *hbar* = *m* = 1. Then we can use a package implementing the Arnoldi method to find the eigenvectors and eigenvalues of *H*.

The script `eigen.jl` looks at the eigen vectors corresponding to the Coulumb potential, *V* ~ 1/*r*.  A script to animate these discretized orbital wavefunctions is `anim.sh`, which can be ran in a terminal as `./anim.sh <res> <N>`, where *N* indicates we will recieve wavefunctions 1 to *N*. If you have anim installed that is.

Here is the short version of what that looks like:

```julia
Δ = laplacian()
V = potential()

H = -0.5Δ + V

decomp = partialschur(H,
                      nev=N,
                      tol=1e-4,
                      which=SR())[1]

E, Ψ = partialeigen(decomp)

```

where the kernel of *V* is:
```julia
coulumb(x, y, z) = -1 / (4π * hypot(x, y, z))
```

## Evolve

