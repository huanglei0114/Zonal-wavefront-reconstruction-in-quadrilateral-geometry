#include "pch.h"
#include "framework.h"
#include "cwfr.h"

// cwfr.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "common.h"


CWFR::CWFR(const MatrixXXd& Sx, const MatrixXXd& Sy, const MatrixXXd& X, const MatrixXXd& Y)
	: m_Sx(Sx)
	, m_Sy(Sy)
	, m_X(X)
	, m_Y(Y)
	, m_rows(Sx.rows())
	, m_cols(Sx.cols())
{
}

CWFR::~CWFR()
{
}

MatrixXXd CWFR::operator()(WFR_METHOD method)
{
	switch (method)
	{
	case WFR_METHOD::HFLI:
		return hfli();
		break;
	case WFR_METHOD::HFLIQ:
		return MatrixXXd();
		break;
	default:
		return MatrixXXd();
		break;
	}
}

MatrixXXd CWFR::hfli()
{
	auto z_size = m_rows * m_cols; // size of the z vector

	/* 0. build the least-squares system */
	/* 0.0 Pre-allocate and reserve the memory space */
	// construct the matrix D and reserve the maximum size
	TripletListd D_trps;
	D_trps.reserve((m_rows - 3) * m_cols + m_rows * (m_cols - 3));

	// construct the std vector g_std for appending the slopes and reserve the maximum size
	std_vecd g_std;
	g_std.reserve((m_rows - 3) * m_cols + m_rows * (m_cols - 3));

	/* 0.1 fill D and g_std */
	hfli_fill_D_g(D_trps, g_std);

	// 1. solve the least-squares system
	// build the sparse matrix D
	SparseMatrixXXd D(D_trps.size() / 2, m_cols * m_rows);
	D.setFromTriplets(D_trps.begin(), D_trps.end());
	D.makeCompressed();

	// map the vecotr g
	VectorMapd g(g_std.data(), g_std.size());

	// solve with QR factorization
	QRSolver qr_solver;
	qr_solver.compute(D);
	if (qr_solver.info() != Eigen::Success) {
		return MatrixXXd();
	}
	
	VectorXd z = qr_solver.solve(g);
	if (qr_solver.info() != Eigen::Success) {
		return MatrixXXd();
	}

	//z.reshaped()
	//// 2. build and return the result
	//MatrixXXd Z = MatrixXXd::Zero(m_rows, m_cols);
	//for (int_t i = 0; i < m_rows; i++) {
	//	for (int_t j = 0; j < m_cols; j++) {
	//		Z(i, j) = z(ID_1D(j, i, m_cols));
	//	}
	//}

	return z.reshaped(m_rows, m_cols);
}

void CWFR::hfli_fill_D_g(TripletListd& D_trps, std_vecd& g_std)
{
	bb_pair is_sx_sy_5th; // determine if 5th order
	bb_pair is_sx_sy_3rd; // determine if 3rd order
	int_t curr_row = 0;

	// start the iterations 
	for (int_t i = 1; i <= m_rows - 2; i++) {
		for (int_t j = 1; j <= m_cols - 2; j++) {
			// validate if 5th,3rd or no equations
			is_sx_sy_5th = is_5th_order_equation(i, j);
			is_sx_sy_3rd = is_3rd_order_equation(i, j);

			// deal with Sx
			if (is_sx_sy_5th.first || is_sx_sy_3rd.first) {
				// push to D_trps
				D_trps.push_back(Tripletd(curr_row, ID_1D(i, j + 1, m_cols),  1));
				D_trps.push_back(Tripletd(curr_row, ID_1D(i, j    , m_cols), -1));
				++curr_row;

				// push to g_std
				if (is_sx_sy_5th.first) g_std.push_back(calculate_5th_order_gx(i, j));
				else g_std.push_back(calculate_3rd_order_gx(i, j));
			}

			// deal with Sy
			if (is_sx_sy_5th.second || is_sx_sy_3rd.second) {
				// push to D_trps
				D_trps.push_back(Tripletd(curr_row, ID_1D(i    , j, m_cols),  1));
				D_trps.push_back(Tripletd(curr_row, ID_1D(i + 1, j, m_cols), -1));
				++curr_row;

				// push_to g_std
				if (is_sx_sy_5th.second) g_std.push_back(calculate_5th_order_gy(i, j));
				else g_std.push_back(calculate_3rd_order_gy(i, j));
			}
		}
	}
}

bb_pair CWFR::is_3rd_order_equation(const int_t& i, const int_t& j)
{
	std::pair<bool, bool> is_valid{ true, true };

	// slope x case
	if (
		!std::isfinite(m_Sx(i, j)) ||
		!std::isfinite(m_Sx(i, j + 1)) ||
		!std::isfinite(m_Sx(i, j - 1))
		)
	{
		is_valid.first = false;
	}

	// slope y case
	if (
		!std::isfinite(m_Sy(i, j)) ||
		!std::isfinite(m_Sy(i + 1, j)) ||
		!std::isfinite(m_Sy(i - 1, j))
		)
	{
		is_valid.second = false;
	}

	return is_valid;
}

bb_pair CWFR::is_5th_order_equation(const int_t& i, const int_t& j)
{
	std::pair<bool, bool> is_valid{ true, true };

	// it should be 3rd-order valid
	is_valid = is_3rd_order_equation(i, j);

	// slope x case
	if (is_valid.first) {
		// it should be valid at the right boundary
		is_valid.first = j + 2 > m_cols - 1 ? false : true;

		// it should avoid the NaN at the j + 2 position
		if (is_valid.first) {
			is_valid.first = std::isfinite(m_Sx(i, j + 2)) ? true : false;
		}
	}

	// slope y case
	if (is_valid.second) {
		// it should be valid at the bottom boundary
		is_valid.second = i + 2 > m_rows - 1 ? false : true;

		// it should avoid the NaN at the i + 2 position
		if (is_valid.second) {
			is_valid.second = std::isfinite(m_Sy(i + 2, j)) ? true : false;
		}
	}

	return is_valid;
}

double CWFR::calculate_3rd_order_gx(const int_t& i, const int_t& j)
{
	return (m_Sx(i, j) + m_Sx(i, j + 1)) * (m_X(i, j + 1) - m_X(i, j)) * 0.5;
}

double CWFR::calculate_5th_order_gx(const int_t& i, const int_t& j)
{
	return (-1 / 13 * m_Sx(i, j - 1) + m_Sx(i, j) + m_Sx(i, j + 1) - 1 / 13 * m_Sx(i, j + 2)) * (m_X(i, j + 1) - m_X(i, j)) * (13 / 24);
}

double CWFR::calculate_3rd_order_gy(const int_t& i, const int_t& j)
{
	return (m_Sy(i, j) + m_Sy(i + 1, j)) * (m_Y(i + 1, j) - m_Y(i, j)) * 0.5;
}

double CWFR::calculate_5th_order_gy(const int_t& i, const int_t& j)
{
	return (-1 / 13 * m_Sy(i - 1, j) + m_Sy(i, j) + m_Sy(i + 1, j) - 1 / 13 * m_Sy(i + 2, j)) * (m_Y(i + 1, j) - m_Y(i, j)) * (13 / 24);
}
