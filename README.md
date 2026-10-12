# Shallow Water Solver — HPSC Project (INFO0939, 2024)

A 2D linear shallow water equations solver written in C, parallelised five ways
for the *High Performance Scientific Computing* course project (INFO0939,
University of Liège). Course and project page:
<https://people.montefiore.uliege.be/geuzaine/INFO0939/project/>.

The solver reads a bathymetry (water depth `h`) on a regular grid, then advances
the free-surface elevation `eta` and velocities `u`, `v` with an explicit
finite-difference scheme on a staggered grid. The elevation field is written as
VTK image files (`.vti` + a `.pvd` time series) that you can open in ParaView.

## Repository layout

```
.
├── Makefile          # builds every version + run shortcuts
├── src/
│   ├── shallow.c          # sequential reference
│   ├── shallow_omp.c      # OpenMP (shared memory)
│   ├── shallow_mpi.c      # MPI (2D Cartesian domain decomposition)
│   ├── shallow_mpi_omp.c  # hybrid MPI + OpenMP
│   └── shallow_gpu.c      # OpenMP target offloading (GPU)
├── params/
│   ├── h_simple.dat             # input bathymetry (binary)
│   ├── param_simple.txt         # default run (800 x 800 grid)
│   ├── param_simple_original.txt
│   └── param_simple_<N>.txt     # weak-scaling inputs for N = 1..64 cores
└── jobs/             # Slurm scripts (strong/weak scaling, Scalasca profiling)
```

## Quick start

Requirements: a C compiler with OpenMP (`gcc`), and an MPI implementation
(OpenMPI) for the MPI versions.

```bash
make shallow          # build the sequential version
make run_sequential   # run it on params/param_simple.txt
```

Output:

```
...
 - grid size: 4000 m x 4000 m (800 x 800 = 640000 grid points)
 - number of time steps: 4000
Done: 2.2 seconds (1163 MUpdates/s)
```

Results land in `run/` (`eta_simple_*.vti`, `eta_simple.pvd`). Open the
`.pvd` file in ParaView to animate the wave. `make clean` removes executables
and `run/`.

## Building and running every version

`make` (or `make all`) builds all five executables at the repository root.

| Version          | Build                  | Run (via make)        | Default setup          |
|------------------|------------------------|-----------------------|------------------------|
| Sequential       | `make shallow`         | `make run_sequential` | 1 core                 |
| OpenMP           | `make shallow_omp`     | `make run_omp`        | 4 threads              |
| MPI              | `make shallow_mpi`     | `make run_mpi`        | 4 ranks                |
| MPI + OpenMP     | `make shallow_mpi_omp` | `make run_mpi_omp`    | 4 ranks × 4 threads    |
| GPU (offloading) | `make shallow_gpu`     | `make run_gpu`        | `OMP_TARGET_OFFLOAD=MANDATORY` |

To run by hand, call any executable with a parameter file. **The bathymetry
file named in the parameter file is opened relative to the current
directory**, so copy it next to where you run:

```bash
mkdir -p run && cp params/h_simple.dat run/ && cd run
OMP_NUM_THREADS=8 ../shallow_omp ../params/param_simple.txt
mpirun -np 8 ../shallow_mpi ../params/param_simple.txt
```

Use a different parameter file with `make run_mpi PARAM_FILE=../params/param_simple_original.txt`.

**GPU note:** `shallow_gpu.c` uses `#pragma omp declare mapper`, which needs
clang or GCC ≥ 15. By default it is built with plain `-fopenmp`, so target
regions fall back to the host CPU (and `make run_gpu` fails if no device is
available). For real offloading, e.g. on an NVIDIA A100 (`sm_80`):

```bash
make shallow_gpu CC=clang GPUFLAGS='-fopenmp --offload-arch=sm_80'
```

**macOS note:** Apple's `gcc` is clang without OpenMP. Only the sequential
version builds out of the box; install Homebrew GCC (e.g. `make CC=gcc-14`)
for the others.

## Parameter file format

One value per line, in this order (`params/param_simple.txt`):

| Line | Value          | Meaning                                                        |
|------|----------------|----------------------------------------------------------------|
| 1    | `5`            | `dx` — grid spacing in x (m)                                   |
| 2    | `5`            | `dy` — grid spacing in y (m)                                   |
| 3    | `0.05`         | `dt` — time step (s)                                           |
| 4    | `200`          | `max_t` — simulated time (s)                                   |
| 5    | `9.81`         | `g` — gravitational acceleration (m/s²)                        |
| 6    | `2e-5`         | `gamma` — dissipation coefficient (1/s)                        |
| 7    | `1`            | source type: `1` = sinusoidal velocity on top boundary, `2` = sinusoidal elevation at domain centre |
| 8    | `40`           | output sampling rate (write every N steps, `0` = no output)     |
| 9    | `h_simple.dat` | input bathymetry file                                           |
| 10–12| `eta_simple` … | output base names for `eta`, `u`, `v`                            |

The grid size follows from the bathymetry extent divided by `dx`, `dy`: halving
the spacing quadruples the number of grid points. The `param_simple_<N>.txt`
files scale the problem with the number of cores `N` (constant work per core)
and disable output, which is what the weak-scaling jobs use.

## Scaling experiments on a Slurm cluster

The scripts in `jobs/` were written for the course cluster. They expect the
repository to be cloned at `~/project_info0939` with the executables built at
its root (edit `PROJECT_DIR` at the top of each script otherwise).

```bash
cd ~/project_info0939
make
sbatch jobs/strong_scaling_mpi.job
```

| Script                       | What it measures                                                    |
|------------------------------|---------------------------------------------------------------------|
| `strong_scaling_openmp.job`  | fixed problem, 1 → all threads, `OMP_PROC_BIND=close/spread`        |
| `strong_scaling_mpi.job`     | fixed problem, 1 → all ranks, ranks packed vs spread                |
| `strong_scaling_hybrid.job`  | fixed problem, 2/4/8 threads per rank, increasing ranks             |
| `weak_scaling_openmp.job`    | problem grows with threads (`param_simple_<N>.txt`)                 |
| `weak_scaling_mpi.job`       | problem grows with ranks                                            |
| `weak_scaling_hybrid.job`    | problem grows with ranks × threads                                  |
| `scalasca_profile.job`       | Scalasca trace of the hybrid version on `param_simple_original.txt` |

Each job prints a timing table to `job_output_<kind>_<jobid>.out` and keeps
per-run outputs in `results_<kind>_<jobid>/`.
