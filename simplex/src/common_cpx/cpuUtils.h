#pragma once

#include <Eigen/Dense>
#include "para.h"

#ifdef USE_DOUBLE
typedef Eigen::VectorXd VECTOR;
typedef Eigen::MatrixXd MATRIX;
typedef Eigen::VectorXi VECTORi;
#else
typedef Eigen::VectorXf VECTOR;
typedef Eigen::MatrixXf MATRIX;
typedef Eigen::VectorXi VECTORi;
#endif