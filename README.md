# Near-Optimal Convergence of Widlund's Method and Numerical Exploration of Different Implementations
This repository is a companion to the paper [ADD LINK WHEN PUBLISHED](https://github.com/LoganOReed/WidlundImplementations). It contains the implementations used for benchmarking, without the tooling needed to report the residuals and iterates at each iteration. These versions of the implementations also *do not* run a final check with the true residual, b-Ax, and trusts the user to ensure proper convergence.

## Usage
The implementations are fully stripped down, so the repository can function both as a way to test the reproducibility of the results and as a collection of implementations which could be easily integrated with other code. Notably, there is no check that the computed residual and the true residual's coincide. This can be easily reintroduced by computing the true residual and checking the result of the function call.

## Authors
The companion paper is written by 
- Andrew Higgins
- Logan Reed
- Daniel Szyld

flexGAL.m is a modified version of [gal_flex.m](https://git.uni-wuppertal.de/anginf/krylov4phs), which was simplified to only consider real vector spaces. 

## License
GPL v3.0
