#ifndef CWFR_H
#define CWFR_H

#ifdef WAVEFRONTRECONSTRUCTION_EXPORTS
#define WAVEFRONTRECONSTRUCTION_API __declspec(dllexport)
#else
#define WAVEFRONTRECONSTRUCTION_API __declspec(dllimport)
#endif

// Eigen API aliases
using int_t = Eigen::Index;
using Tripletd = Eigen::Triplet<double>;
using SparseMatrixXXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using MatrixXXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXd = Eigen::Vector<double, Eigen::Dynamic>;


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

public:
	CWFR() = default;
	virtual ~CWFR();

	// Disable copyping
	CWFR(const CWFR&) = delete;
	CWFR& operator=(const CWFR&) = delete;

	MatrixXXd operator () (
		const MatrixXXd& Sx,/*!< [in] Slopes in x direction*/
		const MatrixXXd& Sy,/*!< [in] Slopes in y direction*/
		const MatrixXXd& X, /*!< [in] x coordinates*/
		const MatrixXXd& Y, /*!< [in] y coordinates*/
		WFR_METHOD method = WFR_METHOD::HFLI /*!< [in] method to be used*/
	);

private:
	//! HFLI method
	/*!
	* Reconstruct the height from the slopes in x and y directions with
	* the High-order Finite-difference-based Least-squares Integration (HFLI)
	* method.
	*					D * z = g
	* \return the reconstructed wavefront Z
	*/
	MatrixXXd hfli(
		const MatrixXXd& Sx,/*!< [in] Slopes in x direction*/
		const MatrixXXd& Sy,/*!< [in] Slopes in y direction*/
		const MatrixXXd& X, /*!< [in] x coordinates*/
		const MatrixXXd& Y  /*!< [in] y coordinates*/
	);

private:
	void preapre_Dzg(
		
	);
};


#endif // !CWFR_H


