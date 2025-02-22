NAME = CurvatureEngine

CC = clang++

CFLAGS = -g -std=c++17 -O3 \
         -mavx2 -mfma -march=native -mtune=native \
         -funroll-loops -fvectorize -ffp-contract=fast \
         -freciprocal-math -ffast-math -fstrict-aliasing \
         -fomit-frame-pointer -flto=full -mprefer-vector-width=256 \
         -I$(INC_DIR)



SRC_DIR = srcs
INC_DIR = includes
OBJ_DIR = build
OUT_VTK_DIR = output


SRC = $(wildcard $(SRC_DIR)/*.cpp)
TOTAL := $(words $(SRC))


OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)


GREEN := \033[0;32m
YELLOW := \033[0;33m
RED := \033[0;31m
NC := \033[0m

.PHONY: all clean fclean re

all: $(NAME)
	@if [ -f $(OBJ_DIR)/.counter ]; then rm $(OBJ_DIR)/.counter; fi
	@echo -e "$(GREEN)Build complete!$(NC)"


$(NAME): $(OBJ)
	@echo -e "$(YELLOW)Linking $@...$(NC)"
	$(CC) $(CFLAGS) -o $@ $^ 


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@if [ ! -f $(OBJ_DIR)/.counter ]; then echo 0 > $(OBJ_DIR)/.counter; fi; \
	cnt=$$(cat $(OBJ_DIR)/.counter); \
	cnt=$$((cnt+1)); \
	echo $$cnt > $(OBJ_DIR)/.counter; \
	perc=$$((100 * cnt / $(TOTAL))); \
	echo -e "$(GREEN)Compiling $< [$$cnt/$(TOTAL) ($$perc%%)]$(NC)"; \
	$(CC) $(CFLAGS) -c $< -o $@ 

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	@echo -e "$(RED)Cleaning up...$(NC)"
	rm -rf $(OBJ_DIR)
	rm -f $(OUT_VTK_DIR)/*.vtk

fclean: clean
	rm -f $(NAME)

re: fclean all
