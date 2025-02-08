NAME = CurvatureEngine
CC = clang++

CFLAGS = -g -std=c++17 -O3 \
		 -mavx2 -mfma -march=native -mtune=native \
		 -funroll-loops -fvectorize -ffp-contract=fast \
		 -freciprocal-math -ffast-math -fstrict-aliasing \
		 -fomit-frame-pointer -flto=full -mprefer-vector-width=256

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(NAME)

fclean: clean

re: fclean all

.PHONY: all clean fclean re

