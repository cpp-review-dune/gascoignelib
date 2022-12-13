#include <gtest/gtest.h>

#include <Mesh/p4estmeshagent.h>
#include <Mesh/p8estmeshagent.h>

using namespace Gascoigne;

TEST(pforest_test, basic_tree_test)
{
  P4estMeshAgent pma("data/square.inp", 0);
  EXPECT_EQ(pma.trees_count(), 1);
  EXPECT_EQ(pma.quad_count(), 1);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 4);
  IndexVector refine_cells;
  refine_cells.push_back(2);
  pma.refine_cells(refine_cells);
  EXPECT_EQ(pma.quad_count(), 7);
}

TEST(pforest_test, multitree_test)
{
  P4estMeshAgent pma("data/square_2.inp", 0);
  EXPECT_EQ(pma.trees_count(), 2);
  EXPECT_EQ(pma.quad_count(), 2);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 8);
}

TEST(pforest_test, threeDtree_test)
{
  P8estMeshAgent pma("data/box.inp", 0);
  EXPECT_EQ(pma.trees_count(), 1);
  EXPECT_EQ(pma.quad_count(), 1);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 8);
}
