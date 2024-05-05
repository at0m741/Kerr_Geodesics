NAME = geodesic_solver
NAME_MPI = geodesic_solver_mpi

CC = gcc
MPI = mpicc

FLAGS = -lm -g -fopenmp -masm=intel -ffast-math -funroll-loops -mavx2
MPICC_FLAGS = -g -masm=intel -ffast-math -funroll-loops -mavx2

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: COMPILER = $(CC)
all: CFLAGS = $(FLAGS)
all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(FLAGS)

$(NAME_MPI): $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_MPI) $(MPICC_FLAGS)

%.o: %.c
	$(COMPILER) -c $< -o $@ $(CFLAGS)

clean:
	rm -f $(OBJ)

fclean: clean
	rm -f $(NAME)

re: fclean all

mpicc: COMPILER = $(MPI)
mpicc: CFLAGS = $(MPICC_FLAGS)
mpicc: $(OBJ)
	$(MPI) $(OBJ) -o $(NAME_MPI) $(MPICC_FLAGS)

.PHONY: all clean fclean re