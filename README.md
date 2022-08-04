# Tree-Structure-Algorithm

Tree structure algorithm (TSA) for finite horizon optimal control problems with Hamilton-Jacobi-Bellman (HJB) equation.

## Contents
### Numerical test scripts
 - `main_git.m` Test script of the tree structure with heat equation generating resolution of HJB and optimal control and trajectory.

### Tree construction
 - `full_tree_pruning.m` Create the tree structure with pruning rule for the given PDE by generating the nodes, length of each level of the tree and the adjacency matrix.
 - `tree_creation_heat.m` Create the tree structure without pruning for the heat equation.

### Auxiliary
- `check.m` Function to deduce if the pruning rule is active.
- `tridiag.m` Function to solve equations with tridiagonal systems to implement the Thomas algorithm.
