include Mk/libs_flags.mk


OBJS=pmm_trajectory.o pmm_trajectory3d.o velocity_search_graph.o three_acc.o
OBJS_TESTS=pmm_trajectory_unittest.o
TARGET=main
TEST_TARGET=unittests
OBJ_DIR=obj

include Mk/recipe.mk
