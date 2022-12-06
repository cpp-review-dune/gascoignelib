#include <gtest/gtest.h>

#include <Mesh/pforestmeshagent.h>

using namespace Gascoigne;

// Demonstrate some basic assertions.
TEST(pforest_test, basic_tree_test)
{
  P4estMeshAgent pma("square.inp", 0);
  EXPECT_EQ(pma.trees_count(), 1);
  EXPECT_EQ(pma.quad_count(), 1);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 4);
  IndexVector refine_cells;
  refine_cells.push_back(2);
  pma.refine_cells(refine_cells);
  EXPECT_EQ(pma.quad_count(), 7);
}

// Demonstrate some basic assertions.
TEST(pforest_test, multitree_test)
{
  P4estMeshAgent pma("square_2.inp", 0);
  EXPECT_EQ(pma.trees_count(), 2);
  EXPECT_EQ(pma.quad_count(), 2);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 8);
}
