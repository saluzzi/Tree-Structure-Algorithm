# Tree-Structure-Algorithm

Tree Structure Algorithm (TSA) for finite horizon optimal control problems with Hamilton-Jacobi-Bellman (HJB) equation. 
Based on the paper: A. Alla, M. Falcone, L. Saluzzi, An efficient DP algorithm on a tree-structure for finite horizon optimal control problems, SIAM J. Sci. Comput., 41 (4), A2384–A2406 , 2019.
(Thanks to Xueying Wang for the help)

## Contents
### Numerical test scripts
 - `main_git.m` Test script of the tree structure driven by the heat equation generating resolution of HJB, optimal control and trajectory.

### Tree construction
 - `full_tree_pruning.m` Create the tree structure with pruning rule for the given dynamics by generating the nodes, length of each level of the tree and the adjacency list.
 - `tree_creation_heat.m` Create the tree structure without pruning for the heat equation.

### Auxiliary
- `check.m` Function to deduce if the pruning rule is active.
- `tridiag.m` Function to solve equations with tridiagonal systems to implement the Thomas algorithm.
