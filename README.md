# Quantum Mechanics Simulations
A collection of scripts to find simulated position space eigen vectors, or wavefunctions, of the discretized time independent and dependent Sch&ouml;dinger Equation, [*ISE*, *DSE*]: 

The time Indpendent (Implicit) Schr&ouml;dinger Equation
```
ISE: H|psi(x)> = E|psi(x)> 
```

The time Dependent (Deterministic) Schr&ouml;dinger Equation
```
DSE: |psi(x, t)> = exp(im*H*t)|psi(x)>
```

## Shoot
The Shooting Method, or "wagging the tail", is an algorithm to numerically solve for the wavefunctions of the *ISE*, for arbitrary potentials in 1-D, including the radial equation for spherically symmetric potentials.

The idea is that the we treat the *ISE* as a boundary value problem, where *psi*(0) = 0 is given, and we want *psi*(R) = 0; we then guess a low value for E, integrate the *ISE* out to *R* = *L*, which is a "shot"; we then see if the *psi*(L) < or > 0; we then step, *E* -> *E* + *dE*, until we cross over 0, and then use bisection to narrow down a solution. When animated it looks like the wavefunction is wagging its tail, hence the popular name. 

A novel approach taken here is that instead of using Euler's method for numerical integration, this simulation implements Numerov's method, a 4th order implicit method suited perfectly for the *ISE*.

I need to figure out how to include `LaTeX` here.

## Eigen
This method "simply" sets up the matrix eigen value problem associated with solving the *ISE* on a 3-D lattice, where each dimension has a resolution of *res* intervals, resulting in a cube in space with *res*^3 lattice points.

This can be flattened into a vector of length *res*^3, with a corresponding Hamiltonian operator, *H* = 0.5*L*^2 + *V*, where *L*^2 is the 3-D Laplacian operator and *hbar* = *m* = 1. Then we can use a package implementing the Arnoldi method to find the eigenvectors an eigenvalues of *H*.

The script `eigen.jl` looks at the eigen vectors corresponding to the Coulumb Potential, *V* ~ 1/*r*^2.  A script to animate these discretized orbital wavefunctions is `anim.sh`, which can be ran in a terminal as `./anim.sh <res> <N>`, where *N* indicates we will recieve wavefunctions 1 to *N*. If you have anim installed that is.



## Evolve

