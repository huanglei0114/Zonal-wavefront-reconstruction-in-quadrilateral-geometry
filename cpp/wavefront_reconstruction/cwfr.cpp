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
		return hfli_calculator(std::bind(&CWFR::hfli_fill_D_g, this, std::placeholders::_1, std::placeholders::_2));
	case WFR_METHOD::HFLIQ:
	default:
		return hfli_calculator(std::bind(&CWFR::hfliq_fill_D_g, this, std::placeholders::_1, std::placeholders::_2));
	}
}

MatrixXXd CWFR::hfli_calculator(std::function<void(TripletListd&, std_vecd&)> hfli_prep)
{
	auto z_size = m_rows * m_cols; // size of the z vector

	/* 0. build the least-squares system */
	/* 0.0 Pre-allocate and reserve the memory space */
	// construct the matrix D and reserve the maximum size
	TripletListd D_trps;
	D_trps.reserve(z_size);

	// construct the std vector g_std for appending the slopes and reserve the maximum size
	std_vecd g_std;
	g_std.reserve(z_size);

	/* 0.1 fill D and g_std */
	hfli_prep(D_trps, g_std);

	/* 1. solve the least - squares system */
	// build the sparse matrix D
	SparseMatrixXXd D(D_trps.size() / 2, z_size);
	D.setFromTriplets(D_trps.begin(), D_trps.end());
	D.makeCompressed();

	// map the vecotr g
	VectorMapd g(g_std.data(), g_std.size());

	// solve with QR factorization
	Solver qr_solver(D);
	VectorXd z = qr_solver.solve(g);
	if (qr_solver.info() != Eigen::Success) {
		return MatrixXXd::Zero(m_rows, m_cols);
	}

	/* 2. Only keep the valid points */
	for (int_t i = 0; i < m_rows; i++) {
		for (int_t j = 0; j < m_cols; j++) {
			if (!isfinite(m_Sx(i, j)) || !isfinite(m_Sy(i, j)))
				z(ID_1D(j, i, m_cols)) = NAN;
		}
	}

	/* 3. Reshape and return the result */
	return z.reshaped<Eigen::RowMajor>(m_rows, m_cols);
}

void CWFR::hfli_fill_D_g(TripletListd& D_trps, std_vecd& g_std)
{
	bool is_5th = false; // determine if 5th order
	bool is_3rd = false; // determine if 3rd order
	int_t curr_row = 0;

	// start the x iterations 
	for (int_t i = 0; i <= m_rows - 1; i++) {
		for (int_t j = 0; j <= m_cols - 2; j++) {

			// validate if 5th,3rd or no equations
			is_5th = is_5th_order_equation_sx(i, j);
			is_3rd = is_3rd_order_equation_sx(i, j);

			// deal with Sx
			if (is_5th || is_3rd) {
				// push to D_trps
				D_trps.push_back(Tripletd(curr_row, ID_1D(j, i, m_cols), -1));
				D_trps.push_back(Tripletd(curr_row, ID_1D(j + 1, i, m_cols), 1));
				++curr_row;

				// push to g_std
				if (is_5th) g_std.push_back(calculate_5th_order_gx(i, j));
				else g_std.push_back(calculate_3rd_order_gx(i, j));
			}
		}
	}

	// start the y iterations 
	for (int_t i = 0; i <= m_rows - 2; i++) {
		for (int_t j = 0; j <= m_cols - 1; j++) {

			// validate if 5th,3rd or no equations
			is_5th = is_5th_order_equation_sy(i, j);
			is_3rd = is_3rd_order_equation_sy(i, j);

			// deal with Sy
			if (is_5th || is_3rd) {
				// push to D_trps
				D_trps.push_back(Tripletd(curr_row, ID_1D(j, i, m_cols), -1));
				D_trps.push_back(Tripletd(curr_row, ID_1D(j, i + 1, m_cols), 1));
				++curr_row;

				// push_to g_std
				if (is_5th) g_std.push_back(calculate_5th_order_gy(i, j));
				else g_std.push_back(calculate_3rd_order_gy(i, j));
			}
		}
	}
}

void CWFR::hfliq_fill_D_g(TripletListd& D_trps, std_vecd& g_std)
{
	bool is_5th = false; // determine if 5th order
	bool is_3rd = false; // determine if 3rd order
	int_t curr_row = 0;

	// start the x iterations 
	for (int_t i = 0; i <= m_rows - 1; i++) {
		for (int_t j = 0; j <= m_cols - 2; j++) {
			// validate if 5th,3rd or no equations
			is_5th = is_5th_order_equation_sx(i, j);
			is_3rd = is_3rd_order_equation_sx(i, j);

			// deal with Sx
			if (is_5th || is_3rd) {
				// push to D_trps
				D_trps.push_back(Tripletd(curr_row, ID_1D(j, i, m_cols), -1));
				D_trps.push_back(Tripletd(curr_row, ID_1D(j + 1, i, m_cols), 1));
				++curr_row;

				// push to g_std
				if (is_5th)
					g_std.push_back(
						calculate_5th_order_gx(i, j, false) + calculate_5th_order_gx(i, j, true)
					);
				else
					g_std.push_back(
						calculate_3rd_order_gx(i, j, false) + calculate_3rd_order_gx(i, j, true)
					);
			}
		}
	}

	// start the y iterations 
	for (int_t i = 0; i <= m_rows - 2; i++) {
		for (int_t j = 0; j <= m_cols - 1; j++) {

			// validate if 5th,3rd or no equations
			is_5th = is_5th_order_equation_sy(i, j);
			is_3rd = is_3rd_order_equation_sy(i, j);

			// deal with Sy
			if (is_5th || is_3rd) {
				// push to D_trps
				D_trps.push_back(Tripletd(curr_row, ID_1D(j, i, m_cols), -1));
				D_trps.push_back(Tripletd(curr_row, ID_1D(j, i + 1, m_cols), 1));
				++curr_row;

				// push_to g_std
				if (is_5th)
					g_std.push_back(
						calculate_5th_order_gy(i, j, false) + calculate_5th_order_gy(i, j, true)
					);
				else
					g_std.push_back(
						calculate_3rd_order_gy(i, j, false) + calculate_3rd_order_gy(i, j, true)
					);
			}
		}
	}

}

bool CWFR::is_3rd_order_equation_sx(const int_t& i, const int_t& j)
{
	bool is_valid = true;

	// slope x case
	if (
		!std::isfinite(m_Sx(i, j)) ||
		!std::isfinite(m_Sx(i, j + 1))
		)
	{
		is_valid = false;
	}

	return is_valid;
}

bool CWFR::is_3rd_order_equation_sy(const int_t& i, const int_t& j)
{
	bool is_valid = true;

	// slope y case
	if (
		!std::isfinite(m_Sy(i, j)) ||
		!std::isfinite(m_Sy(i + 1, j))
		)
	{
		is_valid = false;
	}

	return is_valid;
}

bool CWFR::is_5th_order_equation_sx(const int_t& i, const int_t& j)
{
	bool is_valid = true;

	// it should be 3rd-order valid
	is_valid = is_3rd_order_equation_sx(i, j);

	// slope x case
	if (is_valid) {
		// it should be valid at the right boundary
		is_valid = (j + 2 <= m_cols - 1) && (j - 1 >= 0) ? true : false;

		// it should avoid the NaN at the j + 2 position
		if (is_valid) {
			is_valid = std::isfinite(m_Sx(i, j + 2)) ? true : false;
		}
	}

	return is_valid;
}

bool CWFR::is_5th_order_equation_sy(const int_t& i, const int_t& j)
{
	bool is_valid = true;

	// it should be 3rd-order valid
	is_valid = is_3rd_order_equation_sy(i, j);

	// slope y case
	if (is_valid) {
		// it should be valid at the bottom boundary
		is_valid = (i + 2 <= m_rows - 1) && (i - 1 >= 0) ? true : false;

		// it should avoid the NaN at the i + 2 position
		if (is_valid) {
			is_valid = std::isfinite(m_Sy(i + 2, j)) ? true : false;
		}
	}

	return is_valid;
}

double CWFR::calculate_3rd_order_gx(const int_t& i, const int_t& j, bool is_from_gy)
{
	if (!is_from_gy)
		return (m_Sx(i, j) + m_Sx(i, j + 1)) * (m_X(i, j + 1) - m_X(i, j)) * 0.5;
	else
		return (m_Sy(i, j) + m_Sy(i, j + 1)) * (m_Y(i, j + 1) - m_Y(i, j)) * 0.5;
}

double CWFR::calculate_5th_order_gx(const int_t& i, const int_t& j, bool is_from_gy)
{
	if (!is_from_gy)
		return (-1.0 / 13.0 * m_Sx(i, j - 1) + m_Sx(i, j) + m_Sx(i, j + 1) - 1.0 / 13.0 * m_Sx(i, j + 2)) * (m_X(i, j + 1) - m_X(i, j)) * (13.0 / 24.0);
	else
		return (-1.0 / 13.0 * m_Sy(i, j - 1) + m_Sy(i, j) + m_Sy(i, j + 1) - 1.0 / 13.0 * m_Sy(i, j + 2)) * (m_Y(i, j + 1) - m_Y(i, j)) * (13.0 / 24.0);
}

double CWFR::calculate_3rd_order_gy(const int_t& i, const int_t& j, bool is_from_gx)
{
	if (!is_from_gx)
		return (m_Sy(i, j) + m_Sy(i + 1, j)) * (m_Y(i + 1, j) - m_Y(i, j)) * 0.5;
	else
		return (m_Sx(i, j) + m_Sx(i + 1, j)) * (m_X(i + 1, j) - m_X(i, j)) * 0.5;
}

double CWFR::calculate_5th_order_gy(const int_t& i, const int_t& j, bool is_from_gx)
{
	if (!is_from_gx)
		return (-1.0 / 13.0 * m_Sy(i - 1, j) + m_Sy(i, j) + m_Sy(i + 1, j) - 1.0 / 13.0 * m_Sy(i + 2, j)) * (m_Y(i + 1, j) - m_Y(i, j)) * (13.0 / 24.0);
	else
		return (-1.0 / 13.0 * m_Sx(i - 1, j) + m_Sx(i, j) + m_Sx(i + 1, j) - 1.0 / 13.0 * m_Sx(i + 2, j)) * (m_X(i + 1, j) - m_X(i, j)) * (13.0 / 24.0);
}
