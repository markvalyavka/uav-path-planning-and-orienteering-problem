include Mk/libs_flags.mk


OBJS=quad_state.o timer.o pmm_trajectory.o pmm_trajectory3d.o velocity_search_graph.o three_acc.o env_config.o misc_helpers.o base_heuristic.o trajectory_calculation.o
OBJS_TESTS=pmm_trajectory_unittest.o
TARGET=main
TEST_TARGET=unittests
OBJ_DIR=obj

include Mk/recipe.mk
