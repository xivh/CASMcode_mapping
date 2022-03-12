#include "casm/mapping/hungarian.hh"

#include "casm/misc/CASM_math.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace CASM::mapping;

TEST(HungarianTest, Test1) {
  Eigen::MatrixXd C(3, 3);
  C << 0., 1., 3.,  //
      2., 1., 0.,   //
      4., 0., 2.;   //

  double cost;
  hungarian::Assignment assignment;
  std::tie(cost, assignment) = hungarian::solve(C);

  EXPECT_TRUE(almost_equal(cost, 0.0));
  EXPECT_EQ(assignment, hungarian::Assignment({0, 2, 1}));
}

TEST(HungarianTest, Test2) {
  Eigen::MatrixXd C(3, 3);
  C << 1., 1., 2.,  //
      2., 1., 1.,   //
      4., 0., 2.;   //

  double cost;
  hungarian::Assignment assignment;
  std::tie(cost, assignment) = hungarian::solve(C);

  EXPECT_TRUE(almost_equal(cost, 2.0));
  EXPECT_EQ(assignment, hungarian::Assignment({0, 2, 1}));
}
