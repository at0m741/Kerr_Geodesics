NAME = geodesic_solver
NAME_MPI = geodesic_solver_mpi
NAME_KNL = geodesic_solver_knl
NAME_PROFILING = geodesic_solver_profiling


echo:
	@echo "SRC = $(SRC)"
	@echo "OBJ = $(OBJ)"

CC = gcc
MPI = mpicc
FLAGS = -lm -O3 -gstabs -fopenmp -masm=intel -ffast-math -funroll-loops -mavx2 -fopt-info-vec-optimized -fopt-info-all  -fsanitize=address -fsanitize=undefined
MPICC_FLAGS = -g -masm=intel -ffast-math -funroll-loops -mavx2 -fopt-info-vec-optimized -fopt-info-all
KNL_FLAGS = -lm -g -fopenmp -masm=intel -ffast-math -funroll-loops -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -fopt-info-vec-optimized -fopt-info-all
PROFILING_FLAGS = -pg -g -gstabs -fopenmp -ffast-math -funroll-loops -mavx2 -lm -fopt-info-vec-optimized	-fopt-info-all
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: COMPILER = $(CC)
all: CFLAGS = $(FLAGS)
all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(FLAGS) &> opti_info.txt
	@echo $(NAME) "compiled with" $(FLAGS)
	@echo "To run the program use: ./$(NAME)"

$(NAME_MPI): $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_MPI) $(MPICC_FLAGS) &> opti_info.txt

$(NAME_KNL): $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_KNL) $(KNL_FLAGS) &> opti_info.txt

$(NAME_PROFILING): $(OBJ)
	$(CC) $(OBJ) -o $(NAME_PROFILING) $(PROFILING_FLAGS) &> opti_info.txt

%.o: %.c
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

fclean_all: fclean fclean_mpi fclean_knl fclean_profiling

re: fclean all

re_mpi: fclean mpicc

re_profiling: fclean_profiling profiling

re_knl: fclean_knl knl


re_knl: fclean $(NAME_KNL)

x86: COMPILER = $(CC)
x86: CFLAGS = $(FLAGS)
x86: COMPILER = $(CC)
x86: CFLAGS = $(FLAGS)
x86: $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(FLAGS) &> optimization_info_x86.txt
	@echo "\n"
	@echo "x86 version compiled"
	@echo $(NAME) "compiled with" $(FLAGS)
	@echo "To run the program use: ./$(NAME)"
	@echo "Optimization information saved in optimization_info_x86.txt"

mpicc: COMPILER = $(MPI)
mpicc: CFLAGS = $(MPICC_FLAGS)
mpicc: $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_MPI) $(MPICC_FLAGS)
	@echo "\n"
	@echo "MPI version compiled"
	@echo $(NAME_MPI) "compiled with" $(MPICC_FLAGS)
	@echo "To run the program use: mpirun -np <number_of_processes> ./$(NAME_MPI)"

knl: COMPILER = $(MPI)
knl: CFLAGS = $(KNL_FLAGS)
knl: $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_KNL) $(KNL_FLAGS)
	@echo "\n"
	@echo "KNL version compiled"
	@echo $(NAME_KNL) "compiled with" $(KNL_FLAGS)
	@echo "To run the program use: mpirun -np <number_of_processes> ./$(NAME_KNL)"

profiling: COMPILER = $(CC)
profiling: CFLAGS = $(PROFILING_FLAGS)
profiling: $(OBJ)
	$(CC) $(OBJ) -o $(NAME_PROFILING) $(PROFILING_FLAGS)
	@echo "\n"
	@echo "Profiling version compiled"
	@echo $(NAME_PROFILING) "compiled with" $(PROFILING_FLAGS)
	@echo "To run the program use: ./$(NAME_PROFILING)"

.PHONY: all clean fclean re x86 mpicc re_mpi profiling re_profiling knl re_knl