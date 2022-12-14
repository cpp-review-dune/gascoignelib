#include <gtest/gtest.h>

#include <Mesh/p4estmeshagent2d.h>
#include <Mesh/p4estmeshagent3d.h>

using namespace Gascoigne;

TEST(pforest_test, basic_tree_test)
{
  P4estMeshAgent2d pma("data/square.inp");
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
  P4estMeshAgent2d pma("data/square_2.inp");
  EXPECT_EQ(pma.trees_count(), 2);
  EXPECT_EQ(pma.quad_count(), 2);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 8);
}

TEST(pforest_test, threeDtree_test)
{
  P4estMeshAgent3d pma("data/box.inp");
  EXPECT_EQ(pma.trees_count(), 1);
  EXPECT_EQ(pma.quad_count(), 1);
  pma.global_refine(1);
  EXPECT_EQ(pma.quad_count(), 8);
}
