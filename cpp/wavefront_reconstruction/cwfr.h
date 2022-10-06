// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the WAVEFRONTRECONSTRUCTION_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// WAVEFRONTRECONSTRUCTION_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
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


// This class is exported from the dll
class WAVEFRONTRECONSTRUCTION_API CWFR {
public:
	CWFR(void);
	virtual ~CWFR();
};

//extern WAVEFRONTRECONSTRUCTION_API int nwavefrontreconstruction;
//
//WAVEFRONTRECONSTRUCTION_API int fnwavefrontreconstruction(void);
