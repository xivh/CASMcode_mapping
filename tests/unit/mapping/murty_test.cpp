#include "casm/mapping/murty.hh"

#include "casm/casm_io/container/stream_io.hh"
#include "casm/global/eigen.hh"
#include "casm/mapping/hungarian.hh"
#include "casm/misc/CASM_math.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace CASM::mapping;

// To help make a new test:
void print_murty_assignments(
    std::vector<std::pair<double, murty::Assignment>> const &assignments) {
  Index i = 0;
  for (auto const &pair : assignments) {
    double const &cost = pair.first;
    murty::Assignment const &assignment = pair.second;
    std::cout << "EXPECT_TRUE(almost_equal(assignments[" << i << "], {" << cost
              << ", {";
    for (auto const &j : assignment) {
      std::cout << j << ", ";
    }
    std::cout << "}}));" << std::endl;
    ++i;
  }
}

bool almost_equal(std::pair<double, murty::Assignment> const &x,
                  std::pair<double, murty::Assignment> const &y) {
  double tol = 1e-5;
  if (almost_equal(x.first, y.first, tol) && x.second == y.second) {
    return true;
  }
  return false;
}

TEST(MurtyTest, Test1) {
  // test exhaustive return of solutions

  Eigen::MatrixXd C(3, 3);
  C << 0., 1., 3.,  //
      2., 1., 0.,   //
      4., 0., 2.;   //

  int k_best = 10;
  auto assignments = murty::solve(hungarian::solve, C, k_best);

  EXPECT_EQ(assignments.size(), 6);
  // clang-format off
  EXPECT_TRUE(almost_equal(assignments[0], {0, {0, 2, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[1], {3, {0, 1, 2, }}));
  EXPECT_TRUE(almost_equal(assignments[2], {5, {1, 2, 0, }}));
  EXPECT_TRUE(almost_equal(assignments[3], {5, {2, 0, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[4], {5, {1, 0, 2, }}));
  EXPECT_TRUE(almost_equal(assignments[5], {8, {2, 1, 0, }}));
  // clang-format on
}

TEST(MurtyTest, Test2) {
  // test min_cost & max_cost

  Eigen::MatrixXd C(3, 3);
  C << 0., 1., 3.,  //
      2., 1., 0.,   //
      4., 0., 2.;   //

  int k_best = 6;
  double min_cost = 2.5;
  double max_cost = 6.5;
  auto assignments =
      murty::solve(hungarian::solve, C, k_best, min_cost, max_cost);

  EXPECT_EQ(assignments.size(), 4);
  // clang-format off
  EXPECT_TRUE(almost_equal(assignments[0], {3, {0, 1, 2, }}));
  EXPECT_TRUE(almost_equal(assignments[1], {5, {1, 2, 0, }}));
  EXPECT_TRUE(almost_equal(assignments[2], {5, {2, 0, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[3], {5, {1, 0, 2, }}));
  // clang-format on
}

TEST(MurtyTest, Test3) {
  // test that greater than k_best solutions are returned if there are ties

  Eigen::MatrixXd C(3, 3);
  C << 0., 1., 3.,  //
      2., 1., 0.,   //
      4., 0., 2.;   //

  int k_best = 3;
  auto assignments = murty::solve(hungarian::solve, C, k_best);

  EXPECT_EQ(assignments.size(), 5);
  // clang-format off
  EXPECT_TRUE(almost_equal(assignments[0], {0, {0, 2, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[1], {3, {0, 1, 2, }}));
  EXPECT_TRUE(almost_equal(assignments[2], {5, {1, 2, 0, }}));
  EXPECT_TRUE(almost_equal(assignments[3], {5, {2, 0, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[4], {5, {1, 0, 2, }}));
  // clang-format on
}

TEST(MurtyTest, Test4) {
  // test cost matrix with negative values

  Eigen::MatrixXd C(3, 3);
  C << 0., 1., 3.,    //
      2., 1., 0.,     //
      -1., -5., -3.;  //

  int k_best = 3;
  auto assignments = murty::solve(hungarian::solve, C, k_best);

  EXPECT_EQ(assignments.size(), 5);
  // clang-format off
  EXPECT_TRUE(almost_equal(assignments[0], {-5., {0, 2, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[1], {-2., {0, 1, 2, }}));
  EXPECT_TRUE(almost_equal(assignments[2], {0., {1, 2, 0, }}));
  EXPECT_TRUE(almost_equal(assignments[3], {0., {2, 0, 1, }}));
  EXPECT_TRUE(almost_equal(assignments[4], {0., {1, 0, 2, }}));
  // clang-format on
}
