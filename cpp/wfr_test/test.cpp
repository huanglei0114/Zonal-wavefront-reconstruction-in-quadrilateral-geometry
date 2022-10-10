#include "pch.h"
#include "common.h"
#include "cwfr.h"

TEST(TestCaseName, TestName) {
  EXPECT_EQ(1, 1);
  EXPECT_TRUE(true);

  MatrixXXd Sx, Sy, X, Y;

  CWFR wfr(Sx, Sy, X, Y);
}