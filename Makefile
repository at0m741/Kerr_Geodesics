NAME = geodesic_solver
NAME_MPI = geodesic_solver_mpi
NAME_KNL = geodesic_solver_knl
NAME_PROFILING = geodesic_solver_profiling
NAME_KNL_GCC = geodesic_solver_knl_gcc

CC = gcc
MPI = mpicc
FLAGS = -g -lm -O3 -Wopenmp-simd -mavx2 -gstabs -ftree-loop-optimize \
		-ftree-loop-distribution -fopenmp -masm=intel -ffast-math -march=native -mtune=native \
		-funroll-loops -mavx2 -fopt-info-all -fopt-info-vec-optimized
MPICC_FLAGS = -g -masm=intel -ffast-math -funroll-loops -mavx2 \
			  -fopt-info-vec-optimized -fopt-info-all -fopenmp -lm -DUSE_MPI -fopenmp 
KNL_FLAGS = -lm -g -fopenmp -masm=intel -ffast-math -funroll-loops \
			-mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl \
			-fopt-info-vec-optimized -fopt-info-all
PROFILING_FLAGS = -pg -g -gstabs -fopenmp -ffast-math -funroll-loops \
				  -mavx2 -lm -fopt-info-vec-optimized -fopt-info-all \
				  -fopt-info-vec-optimized -fopt-info-all
KNL_GCC_FLAGS = -lm -O3 -g -fopenmp -masm=intel -ffast-math -funroll-loops \
			    -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl \
				-fopt-info-vec-optimized -fopt-info-all
SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c,objs/%.o,$(SRC))

$(shell export HDF5_DIR=/nfs/homes/ltouzali/.local/include/hdf5)
$(OBJ): | objs

objs:
	mkdir -p objs

all: COMPILER = $(CC)
all: CFLAGS = $(FLAGS)
all: $(NAME)

$(NAME): objs $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(FLAGS)
	@echo $(NAME) "compiled with" $(FLAGS)
	@echo "To run the program use: ./$(NAME)"

$(NAME_MPI): objs $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_MPI) $(MPICC_FLAGS)

$(NAME_KNL): objs $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_KNL) $(KNL_FLAGS)

$(NAME_PROFILING): objs $(OBJ)
	$(CC) $(OBJ) -o $(NAME_PROFILING) $(PROFILING_FLAGS)

$(NAME_KNL_GCC): objs $(OBJ)
	$(CC) $(OBJ) -o $(NAME_KNL_GCC) $(KNL_GCC_FLAGS)

objs/%.o: src/%.c
	$(COMPILER) -c $< -o $@ $(CFLAGS)

clean:
	rm -f $(OBJ)

fclean: clean
	rm -f $(NAME) $(wildcard *.vtk)

fclean_mpi: clean
	rm -f $(NAME_MPI) $(wildcard *.vtk)

fclean_knl: clean
	rm -f $(NAME_KNL) $(wildcard *.vtk)

fclean_profiling: clean
	rm -f $(NAME_PROFILING) $(wildcard *.vtk) gmon.out
	rm -f $(wildcard *.txt)
	@echo "Profiling files removed"

fclean_knl_gcc: clean
	rm -f $(NAME_KNL_GCC) $(wildcard *.vtk)

fclean_all: fclean fclean_mpi fclean_knl fclean_profiling fclean_knl_gcc
	@echo "All files removed"

re: fclean all

re_mpi: fclean mpicc

re_profiling: fclean_profiling profiling

re_knl: fclean_knl knl

re_knl_gcc: fclean_knl_gcc $(NAME_KNL_GCC)

x86: COMPILER = $(CC)
x86: CFLAGS = $(FLAGS)
x86: objs $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(FLAGS)
	@echo "\n"
	@echo "x86 version compiled"
	@echo $(NAME) "compiled with" $(FLAGS)
	@echo "To run the program use: ./$(NAME)"
	@echo "Optimization information saved in optimization_info_x86.txt"

mpicc: COMPILER = $(MPI)
mpicc: CFLAGS = $(MPICC_FLAGS) 
mpicc: objs $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_MPI) $(MPICC_FLAGS)
	@echo "\n"  
	@echo "MPI version compiled"
	@echo $(NAME_MPI) "compiled with" $(MPICC_FLAGS)
	@echo "To run the program use: mpirun -np <number_of_processes> ./$(NAME_MPI)"

knl: COMPILER = $(MPI)  
knl: CFLAGS = $(KNL_FLAGS)
knl: objs $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_KNL) $(KNL_FLAGS)
	@echo "\n"
	@echo "KNL version compiled" 
	@echo $(NAME_KNL) "compiled with" $(KNL_FLAGS)
	@echo "To run the program use: mpirun -np <number_of_processes> ./$(NAME_KNL)"

profiling: COMPILER = $(CC)
profiling: CFLAGS = $(PROFILING_FLAGS)
profiling: objs $(OBJ)  
	$(CC) $(OBJ) -o $(NAME_PROFILING) $(PROFILING_FLAGS)
	@echo "\n"
	@echo "Profiling version compiled"
	@echo $(NAME_PROFILING) "compiled with" $(PROFILING_FLAGS) 
	@echo "To run the program use: ./$(NAME_PROFILING)"

knl_gcc: COMPILER = $(CC)
knl_gcc: CFLAGS = $(KNL_GCC_FLAGS)
knl_gcc: objs $(OBJ)
	$(CC) $(OBJ) -o $(NAME_KNL_GCC) $(KNL_GCC_FLAGS)  
	@echo "\n"
	@echo "KNL version compiled with GCC"
	@echo $(NAME_KNL_GCC) "compiled with" $(KNL_GCC_FLAGS)
	@echo "To run the program use: ./$(NAME_KNL_GCC)"  

echo:
	@echo "\n COMPILER OPTIONS :\n"
	@echo "-compile the program with the default compiler use: make all" 
	@echo "-compile the program with the MPI compiler use: make mpicc"
	@echo "-compile the program with the KNL compiler use: make knl"  
	@echo "-compile the program with the profiling flags use: make profiling"
	@echo "-compile the program for KNL or Intel Scalable and GCC use: make knl_gcc"
	@echo "-recompile the program use: make re"
	@echo "-recompile the program with the MPI compiler use: make re_mpi"
	@echo "-recompile the program with the KNL compiler use: make re_knl"
	@echo "-recompile the program with the profiling flags use: make re_profiling"  
	@echo "-recompile the program for KNL or Intel Scalable and GCC use: make re_knl_gcc"
	@echo "-recompile the program with the x86 flags use: make x86"  
	@echo "-remove all files use: make fclean_all"

.PHONY: all clean fclean re mpicc knl profiling knl_gcc re_mpi re_knl re_profiling re_knl_gcc x86 fclean_mpi fclean_knl fclean_profiling fclean_knl_gcc fclean_all echo
