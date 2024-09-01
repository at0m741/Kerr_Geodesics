NAME = geodesic_solver
NAME_PROFILING = geodesic_solver_profiling
NAME_KNL_GCC = geodesic_solver_knl_gcc

CC = gcc

AVX2_FLAGS = -g -lm -O3 -Wopenmp-simd -mavx2 -gstabs -ftree-loop-optimize \
		-ftree-loop-distribution -fopenmp -masm=intel -ffast-math -march=native -mtune=native \
		-funroll-loops -mavx2 -flto -falign-functions=32 -fsanitize=address -DAVX2

AVX512_FLAGS = -g -lm -O3 -Wopenmp-simd -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl \
		-gstabs -ftree-loop-optimize -ftree-loop-distribution -fopenmp -masm=intel -ffast-math \
		-march=native -mtune=native -funroll-loops -flto -falign-functions=32 -fsanitize=address -DAVX512F

PROFILING_FLAGS = -pg -g -gstabs -fopenmp -ffast-math -funroll-loops \
				  -mavx2 -lm -fopt-info-vec-optimized -fopt-info-all \
				  -fopt-info-vec-optimized -fopt-info-all -flto

KNL_GCC_FLAGS = -lm -O3 -g -fopenmp -masm=intel -ffast-math -funroll-loops \
			    -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl \
				-fopt-info-vec-optimized -fopt-info-all -flto

ifeq ($(DEBUG), true)
	AVX2_FLAGS += -D DEBUG
	AVX512_FLAGS += -D DEBUG
endif

ifeq ($(AVX2), true)
	AVX2_FLAGS += -D AVX2
endif

ifeq ($(AVX2), true)
	AVX512_FLAGS += -D AVX512
endif

SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c,objs/%.o,$(SRC))

$(shell export HDF5_DIR=/nfs/homes/ltouzali/.local/include/hdf5)
$(OBJ): | objs

objs:
	mkdir -p objs

all: avx2

avx2: COMPILER = $(CC)
avx2: CFLAGS = $(AVX2_FLAGS)
avx2: $(NAME)_avx2

$(NAME)_avx2: objs $(OBJ)
	$(COMPILER) $(OBJ) -o $(NAME) $(CFLAGS)
	@echo $(NAME) "compiled with AVX2 and" $(CFLAGS)
	@echo "To run the program use: ./$(NAME)"

avx512: COMPILER = $(CC)
avx512: CFLAGS = $(AVX512_FLAGS)
avx512: $(NAME)_avx512

$(NAME)_avx512: objs $(OBJ)
	$(COMPILER) $(OBJ) -o $(NAME) $(CFLAGS)
	@echo $(NAME) "compiled with AVX512 and" $(CFLAGS)
	@echo "To run the program use: ./$(NAME)"

profiling: COMPILER = $(CC)
profiling: CFLAGS = $(PROFILING_FLAGS)
profiling: objs $(OBJ)  
	$(COMPILER) $(OBJ) -o $(NAME_PROFILING) $(PROFILING_FLAGS)
	@echo "\n"
	@echo "Profiling version compiled"
	@echo $(NAME_PROFILING) "compiled with" $(PROFILING_FLAGS) 
	@echo "To run the program use: ./$(NAME_PROFILING)"

knl_gcc: COMPILER = $(CC)
knl_gcc: CFLAGS = $(KNL_GCC_FLAGS)
knl_gcc: objs $(OBJ)
	$(COMPILER) $(OBJ) -o $(NAME_KNL_GCC) $(KNL_GCC_FLAGS)  
	@echo "\n"
	@echo "KNL version compiled with GCC"
	@echo $(NAME_KNL_GCC) "compiled with" $(KNL_GCC_FLAGS)
	@echo "To run the program use: ./$(NAME_KNL_GCC)"  

objs/%.o: src/%.c
	$(COMPILER) -c $< -o $@ $(CFLAGS)

clean:
	rm -f $(OBJ)

fclean: clean
	rm -f $(NAME) $(wildcard *.vtk)

fclean_profiling: clean
	rm -f $(NAME_PROFILING) $(wildcard *.vtk) gmon.out
	rm -f $(wildcard *.txt)
	@echo "Profiling files removed"

fclean_knl_gcc: clean
	rm -f $(NAME_KNL_GCC) $(wildcard *.vtk)

fclean_all: fclean fclean_mpi fclean_knl fclean_profiling fclean_knl_gcc
	@echo "All files removed"

re: fclean avx2

re_profiling: fclean_profiling profiling

re_knl_gcc: fclean_knl_gcc $(NAME_KNL_GCC)

echo:
	@echo "\n COMPILER OPTIONS :\n"
	@echo "-compile the program with AVX2 instructions: make avx2"
	@echo "-compile the program with AVX512 instructions: make avx512"
	@echo "-compile the program with the profiling flags use: make profiling"
	@echo "-compile the program for KNL or Intel Scalable and GCC use: make knl_gcc"
	@echo "-recompile the program use: make re"
	@echo "-recompile the program avec les flags de profiling use: make re_profiling"  
	@echo "-recompile the program for KNL or Intel Scalable and GCC use: make re_knl_gcc"
	@echo "-remove all files use: make fclean_all"

.PHONY: all avx2 avx512 clean fclean re profiling knl_gcc fclean_profiling fclean_knl_gcc fclean_all echo
