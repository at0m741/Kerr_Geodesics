NAME = geodesics_solver

CC = gcc

FLAGS = -lm -g -fopenmp -masm=intel -ffast-math -funroll-loops -mavx2

SRC = $(wildcard *.c)

OBJ = $(SRC:.c=.o)

all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(FLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(FLAGS)

clean:
	rm -f $(OBJ)

fclean: clean
	rm -f $(NAME)	

re: fclean all

.PHONY: all clean fclean re