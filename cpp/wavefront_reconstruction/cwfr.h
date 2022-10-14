#ifndef CWFR_H
#define CWFR_H

#include "common.h"

//! This is the class the reconstruct the wavefront shape from gradient data
/*!
* Reference:
	[1] Guanghui Li, Yanqiu Li, Ke Liu, Xu Ma, and Hai Wang, "Improving
	wavefront reconstruction accuracy by using integration equations with
	higher-order truncation errors in the Southwell geometry," J. Opt. Soc.
	Am. A 30, 1448-1459 (2013)
	[2] Lei Huang, Junpeng Xue, Bo Gao, Chao Zuo, and Mourad Idir,"Zonal
	wavefront reconstruction in quadrilateral geometry for phase measuring
	deflectometry," Appl. Opt. 56, 5139-5144 (2017)
*/
class WAVEFRONTRECONSTRUCTION_API CWFR {
public:
	enum class WFR_METHOD {
		HFLI,
		HFLIQ,
	};

private:
	MatrixXXd m_Sx;
	MatrixXXd m_Sy;
	MatrixXXd m_X;
	MatrixXXd m_Y;
	int_t m_rows;
	int_t m_cols;

public:
	CWFR(
		const MatrixXXd& Sx,/*!< [in] Slopes in x direction*/
		const MatrixXXd& Sy,/*!< [in] Slopes in y direction*/
		const MatrixXXd& X, /*!< [in] x coordinates*/
		const MatrixXXd& Y  /*!< [in] y coordinates*/
	);
	virtual ~CWFR();

	// Disable default constructor and copyping
	CWFR() = delete;
	CWFR(const CWFR&) = delete;
	CWFR& operator=(const CWFR&) = delete;

	MatrixXXd operator () (
		WFR_METHOD method = WFR_METHOD::HFLI /*!< [in] method to be used*/
		);

private:
	//! HFLI method
	/*!
	* Reconstruct the height from the slopes in x and y directions with
	* the High-order Finite-difference-based Least-squares Integration for
	* Quadrileteral (HFLIQ) method.
	*					D * z = g
	* \return the reconstructed wavefront Z
	*/
	MatrixXXd hfli_calculator(std::function<void (TripletListd&, std_vecd&)>hfli_prep);

private:
	//! Fill the matrix D and the rhs vector g for hfli
	void hfli_fill_D_g(
		TripletListd& D_trps, /*!< [out] the filled matrix D*/
		std_vecd& g_std /*!< [out] the filled vector g_std*/
	);

	//! Fill the matrix D and the rhs vector g for hfliq
	void hfliq_fill_D_g(
		TripletListd& D_trps, /*!< [out] the filled matrix D*/
		std_vecd& g_std /*!< [out] the filled vector g_std*/
	);

	//! Determine if the current position is valid for a 3rd-order equation
	/*
	* Determine if the current [i, j] can be used to insert a 3rd-order equation to
	* the matrix D and an element into the rhs vector g. The rules are:
	* 1) For Sx,
	*			x, NaN, x
	* 2) For Sy,
	*			 x
	*			NaN
	*			 x
	* where "x" and "NaN" are the invalid positions.
	* \return a pair of (bool, bool) indicating whether the (x, y) has the 3rd-order equation.
	*/
	bool is_3rd_order_equation_sx(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j  /*!< [in] the id in x-axis*/
	);
	bool is_3rd_order_equation_sy(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j  /*!< [in] the id in x-axis*/
	);

	//! Determine if the current position can write a 5th-order equation
	/*
	* Determine if the current [i, j] can be used to write a 5th-order equation
	* in the matrix D. The rules are:
	* 1) For Sx,
	*			   |
	*			x, x, NaN, x, so the right boundary needs two pixels at least
	* 2) For Sy,
	*			 x
	*			 x <--
	*			NaN
	*			 x
	*	so the bottom boundary needs two pixels at least
	* where "x" and "NaN" are the invalid positions for a 5th-order equation.
	* This means 5th-order = 3rd-order + boundary care + NaN
	* \return a pair of (bool, bool) indicating whether the (x, y) has the 5th-order equation.
	*/
	bool is_5th_order_equation_sx(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j  /*!< [in] the id in x-axis*/
	);
	bool is_5th_order_equation_sy(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j  /*!< [in] the id in x-axis*/
	);

	//! Calcualte the 3rd order for the gx
	/*
	* \return the gx for the j, j+1 segment
	*/
	double calculate_3rd_order_gx(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j, /*!< [in] the id in x-axis*/
		bool is_from_gy = false /*!< [in] for hfliq, if gy contributes*/
	);

	//! Calcualte the 5th order for the gx
	/*
	* \return the gx for the j, j+1 segment
	*/
	double calculate_5th_order_gx(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j, /*!< [in] the id in x-axis*/
		bool is_from_gy = false /*!< [in] for hfliq, if gy contributes*/
	);

	//! Calcualte the 3rd order for the gy
	/*
	* \return the gy for the i, i+1 segment
	*/
	double calculate_3rd_order_gy(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j, /*!< [in] the id in x-axis*/
		bool is_from_gx = false /*!< [in] for hfliq, if gx contributes*/
	);

	//! Calcualte the 5th order for the gy
	/*
	* \return the gy for the i, i+1 segment
	*/
	double calculate_5th_order_gy(
		const int_t& i, /*!< [in] the id in y-axis*/
		const int_t& j, /*!< [in] the id in x-axis*/
		bool is_from_gx = false /*!< [in] for hfliq, if gx contributes*/
	);
};



#endif // !CWFR_H
