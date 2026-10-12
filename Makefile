# Compiler settings
CC=gcc
MPICC = mpicc
CFLAGS= -O3 -Wall -Wextra
LDFLAGS= -lm
OMPFLAGS= -fopenmp
# Without offload support, target regions fall back to the host. For real GPU
# offloading (needs clang or GCC >= 15 for declare mapper), override e.g.:
#   make shallow_gpu CC=clang GPUFLAGS='-fopenmp --offload-arch=sm_80'
GPUFLAGS ?= $(OMPFLAGS)

# MPI settings (add these lines)
MPI_INC = $(shell mpicc --showme:compile)
MPI_LIB = $(shell mpicc --showme:link)

# Directories
SRC_DIR = src
PARAM_DIR = params
RUN_DIR = run

# Targets executables
SEQUENTIAL_TGT = shallow
MPI_TGT = shallow_mpi
OMP_TGT = shallow_omp
MPI_OMP_TGT = shallow_mpi_omp
GPU_TGT = shallow_gpu

# Sources files
SEQUENTIAL_SRC = $(SRC_DIR)/shallow.c
MPI_SRC = $(SRC_DIR)/shallow_mpi.c
OMP_SRC = $(SRC_DIR)/shallow_omp.c
MPI_OMP_SRC = $(SRC_DIR)/shallow_mpi_omp.c
GPU_SRC = $(SRC_DIR)/shallow_gpu.c

# Run settings (paths relative to RUN_DIR)
PARAM_FILE = ../$(PARAM_DIR)/param_simple.txt

# Default target
all: $(SEQUENTIAL_TGT) $(OMP_TGT) $(MPI_TGT) $(MPI_OMP_TGT) $(GPU_TGT)

# Rules for building executables
$(SEQUENTIAL_TGT): $(SEQUENTIAL_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(MPI_TGT): $(MPI_SRC)
	$(MPICC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(OMP_TGT): $(OMP_SRC)
	$(CC) $(CFLAGS) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

$(MPI_OMP_TGT): $(MPI_OMP_SRC)
	$(MPICC) $(CFLAGS) $(MPI_INC) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

$(GPU_TGT): $(GPU_SRC)
	$(CC) $(CFLAGS) $(GPUFLAGS) $^ -o $@ $(LDFLAGS)

# Runs happen in RUN_DIR: the parameter file names the bathymetry file
# (h_simple.dat) relative to the working directory, and outputs land there too.
$(RUN_DIR):
	mkdir -p $(RUN_DIR)
	cp $(PARAM_DIR)/h_simple.dat $(RUN_DIR)/

# Rules for running executables
run_sequential: $(SEQUENTIAL_TGT) | $(RUN_DIR)
	cd $(RUN_DIR) && ../$(SEQUENTIAL_TGT) $(PARAM_FILE)

run_mpi: $(MPI_TGT) | $(RUN_DIR)
	cd $(RUN_DIR) && mpirun -np 4 ../$(MPI_TGT) $(PARAM_FILE)

run_omp: $(OMP_TGT) | $(RUN_DIR)
	cd $(RUN_DIR) && OMP_NUM_THREADS=4 ../$(OMP_TGT) $(PARAM_FILE)

run_mpi_omp: $(MPI_OMP_TGT) | $(RUN_DIR)
	cd $(RUN_DIR) && OMP_NUM_THREADS=4 mpirun -np 4 ../$(MPI_OMP_TGT) $(PARAM_FILE)

run_gpu: $(GPU_TGT) | $(RUN_DIR)
	cd $(RUN_DIR) && OMP_TARGET_OFFLOAD=MANDATORY ../$(GPU_TGT) $(PARAM_FILE)

# Clean up build artifacts and run outputs
clean:
	rm -f $(SEQUENTIAL_TGT) $(MPI_TGT) $(OMP_TGT) $(MPI_OMP_TGT) $(GPU_TGT)
	rm -rf $(RUN_DIR)

.PHONY: all clean run_sequential run_mpi run_omp run_mpi_omp run_gpu
