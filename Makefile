include Mk/libs_flags.mk

OBJS=pmm_trajectory.o pmm_trajectory3d.o velocity_search_graph.o three_acc.o main.o
TARGET=main
OBJ_DIR=obj

include Mk/recipe.mk
