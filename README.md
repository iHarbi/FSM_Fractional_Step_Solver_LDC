# Lid-driven Cavity Benchmark

## Overview

The lid-driven cavity benchmark is a widely recognized fundamental test for Computational Fluid Dynamics (CFD) algorithms. Originally introduced by Odus R. Burggraf in 1966, this scenario involves a fluid contained within a square cavity, where Dirichlet boundary conditions of velocity are applied to all four sides of the cavity.

## Problem Description

In this benchmark:
- The cavity has equal height and length (\(L = H\)).
- The top wall of the cavity moves with a constant velocity of \(|\mathbf{U}| = 1\).
- The remaining three walls are stationary.

The benchmark is used to evaluate the laminar steady-state flow that develops within the cavity, with characteristics that can be compared against a substantial body of established data.

## Historical Context

The lid-driven cavity problem was first proposed by Odus R. Burggraf in 1966 to provide an analytical solution for flow within a cavity. Over time, this scenario has become a standard validation test for CFD algorithms. The most referenced dataset for this benchmark comes from the work of Ghia et al. in 1982, which offers detailed insights into velocity distributions, streamline patterns, and vorticity contours across varying Reynolds numbers.

## Validation and Comparison

The benchmark solution converges to a laminar steady state, and the resulting flow characteristics can be compared to the canonical data documented by Ghia et al. This data serves as a reference for validating CFD simulations, ensuring accuracy in the prediction of velocity profiles and flow patterns within the cavity.

## References

1. Burggraf, Odus R. "Analytical and Numerical Studies of the Structure of Steady Separated Flows." Journal of Fluid Mechanics, vol. 24, no. 1, 1966, pp. 113-151.
2. Ghia, U., Ghia, K. N., and Shin, C. T. "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method." Journal of Computational Physics, vol. 48, no. 3, 1982, pp. 387-411.
