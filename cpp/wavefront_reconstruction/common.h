#ifndef COMMON_H
#define COMMON_H

#ifdef WAVEFRONTRECONSTRUCTION_EXPORTS
#define WAVEFRONTRECONSTRUCTION_API __declspec(dllexport)
#else
#define WAVEFRONTRECONSTRUCTION_API __declspec(dllimport)
#endif

// Aliases
using int_t = Eigen::Index;
using Tripletd = Eigen::Triplet<double>;
using TripletListd = std::vector<Tripletd>;
using SparseMatrixXXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using MatrixXXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;
using std_vecd = std::vector<double>;
using bb_pair = std::pair<bool, bool>;
#define REP(i, a, b) for(int_t i = int_t(a); i <= int_t(b); i++);

#endif // !COMMON_H


