## An SDP-based Branch-and-Cut Algorithm for Biclustering

**BICL-SDP** is an exact algorithm, based on the branch-and-cut technique, for solving the biclustering problem through the $k$-densest-disjoint biclique criterion described in the paper ["An SDP-based Branch-and-Cut Algorithm for Biclustering"](https://arxiv.org/abs/2403.11351). This repository contains the C++ source code, the MATLAB scripts, and the datasets used for the experiments.

## Installation
**BICL-SDP** calls the semidefinite programming solver [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/) by using the [MATLAB Engine API](https://www.mathworks.com/help/matlab/calling-matlab-engine-from-cpp-programs.html) for C++. It requires the MATLAB engine library *libMatlabEngine* and the Matlab Data Array library *libMatlabDataArray*. **BICL-SDP** calls the linear programming solver [Gurobi](https://www.gurobi.com/) and uses [Armadillo](http://arma.sourceforge.net/) to handle matrices and linear algebra operations efficiently.


Ubuntu and Debian instructions:

1) Install MATLAB (>= 2021b)

2) Install Gurobi (>= 10.0.2)

3) Install CMake, OpenBLAS, LAPACK and Armadillo:
 ```
sudo apt-get update
sudo apt-get install cmake libopenblas-dev liblapack-dev libarmadillo-dev
```
4) Open the makefile `biclustering_cpp/Makefile` 
	- Set the variable `matlab_path` with your MATLAB folder.

5) Compile the code:

```
cd biclustering_cpp/
make
```

4) Download SDPNAL+, move the folder `biclustering_matlab` containing the MATLAB source code of **BICL-SDP** in the SDPNAL+ main directory and set the parameter `SDP_SOLVER_FOLDER` of the configuration file accordingly. This folder and its subfolders will be automatically added to the MATLAB search path when **BICL-SDP** starts.

This code has been tested under Ubuntu 22.04 LTS with MATLAB R2021b, Gurobi 10.0.2 and Armadillo 12.6.

## Configuration
Various parameters used in **BICL-SDP** can be modified in the configuration file `biclustering_cpp/config.txt`:

- `BRANCH_AND_BOUND_TOL` - optimality tolerance of the exact algorithm
- `BRANCH_AND_BOUND_PARALLEL` -  thread pool size: single thread (1), multi-thread (> 1)
- `BRANCH_AND_BOUND_MAX_NODES` - maximum number of nodes
- `BRANCH_AND_BOUND_VISITING_STRATEGY` - best first (0),  depth first (1), breadth first (2)
- `MATLAB_SESSION_THREADS_ROOT` - number of threads for the MATLAB session at the root noee
- `MATLAB_SESSION_THREADS_CHILD` - number of threads for the MATLAB session for child nodes
- `SDP_SOLVER_FOLDER` - full path of SDPNAL+ folder
- `SDP_SOLVER_TOL` - accuracy of SDPNAL+ in the relative KKT residual
- `SDP_SOLVER_VERBOSE` - do not display log (0), display log (1)
- `CP_MAX_ITER` - maximum number of cutting-plane iterations
- `CP_TOL` - tolerance between two consecutive cutting-plane iterations
- `CP_MAX_INEQ` - maximum number of valid inequalities to separate
- `CP_PERC_INEQ` - fraction of the most violated inequalities to add
- `CP_EPS_INEQ` - tolerance for checking the violation of the inequalities
- `CP_EPS_ACTIVE` - tolerance for detecting active inequalities
- `GUROBI_FOLDER` - Gurobi solver path
- `GUROBI_VERBOSE` - do not display log (0), display log (1)

## Usage
```
cd biclustering_cpp/
./bb <W_PATH> <K> <LOG_PATH> <RESULT_PATH>
```
- `W_PATH` - full path of the data matrix
- `K` - number of biclusters
- `LOG_PATH` - path of the log file
- `RESULT_PATH` - path of the optimal bicluster assignment matrices

File `W_PATH` contains the weights `w_ij` and the must include an header line with the number of rows `n` and columns `m`:

```
n m
w_11 w_12 ... w_1m
w_21 w_22 ... w_2m
...
...
w_n1 w_n2 ... w_nm
```

## Log

The log file reports the progress of the algorithm:

- `N` - number of rows at the current node
- `M` - number of columns at the current node
- `ID_PAR` - id of the parent node
- `ID` - id of the current node
- `UB_PAR` - upper bound of the parent node
- `UB` - upper bound of the current node
- `TIME (s)` - running time in seconds of the current node
- `CP_ITER` - number of cutting-plane iterations
- `CP_FLAG` - termination flag of the cutting-plane procedure
    - `-2` - maximum number of iterations
    - `-1` - SDP not solved or partially solved successfully
    -  ` 0` - no violated inequalities
    -  ` 1` - node must be pruned
    -  ` 2` - upper bound greater than the previous one
    -  ` 3` - upper bound decrease is not sufficiently large
- `CP_INEQ` - number of inequalities added in the last cutting-plane iteration
- `LB` - current lower bound
- `BEST_LB` - global lower bound
- `SET` - vertex set selection for branching
    -  `U` - branch on the vertices in U
    -  `V` - branch on the vertices in V
    -  `-1` - branching is not needed
- `I J` - indices of branching decision
- `NODE_GAP` - gap at the current node
- `GAP` - overall gap 
- `OPEN` - number of open nodes
