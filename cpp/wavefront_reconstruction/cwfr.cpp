// cwfr.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "cwfr.h"


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
	SparseMatrixXXd D;
	D.reserve(VectorXi::Constant(m_rows, 2));

	// construct the std vector g_std for appending the slopes and reserve the maximum size
	std::vector<double> g_std;
	g_std.reserve((m_rows - 3) * m_cols + m_rows * (m_cols - 3));

	/* 0.1 fill D and g_std */
	hfli_fill_D_g(D, g_std);

	// 1. solve the least-squares system

	// 2. build and return the result
	return MatrixXXd();
}

void CWFR::hfli_fill_D_g(SparseMatrixXXd& D, std::vector<double>& g_std)
{

}

std::pair<bool, bool> CWFR::is_3rd_order_equation(const int_t& i, const int_t& j)
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

std::pair<bool, bool> CWFR::is_5th_order_equation(const int_t& i, const int_t& j)
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
