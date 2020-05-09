# Quantum Mechanics Simulations
A collection of scripts to find simulated position space eigen vectors, or wavefunctions, of the discretized time independent and dependent Sch&ouml;dinger Equation, [ISE, DSE]: 

```
ISE: H*psi(x) = E*psi(x) 

DSE: psi(x, t) = exp(im*H*t)*psi(x)
```

## Shoot
The Shooting Method, or "wagging the tail", is an algorithm to numerically solve for the wavefunctions of the ISE, for arbitrary potentials in 1-D, including the radial equation for spherically symmetric potentials.

The idea is that the we treat the ISE as a boundary value problem, where psi(0) = 0 is given, and we want psi(R) = 0; we then guess a low value for E, integrate the ISE out to R = L, which is a "shot"; we then see if the psi(L) < or > 0; we then step through E, E -> E + dE, until we cross over 0, and then use bisection to narrow down a solution. When animated it looks like the wavefunction is wagging its tail, hence the popular name. 

An added twist is that instead of using Euler's method for numerical integration, this simulation utilizes Numerov's method, a 4th order implicit method suited for equations of the form we are looking at here.

I need to figure out how to include LaTeX at this point.

## Eigen
This method "simply" sets up the matrix eigen value problem associated with solving the ISE on a 3-D lattice, where each dimension has been divided into *d* intervals, resulting in a cube in space with *d*^3 lattice points.  This can be flattened into a matrix of  



## Evolve

## Paths

