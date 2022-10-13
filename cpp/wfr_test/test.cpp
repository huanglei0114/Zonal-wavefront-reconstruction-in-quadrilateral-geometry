#include "pch.h"
#include "common.h"
#include "cwfr.h"
#include "matrix_io.h"

TEST(MatrixIOTest, ReadTheMatrix) {
	const char* file_name = "../../data/X.bin";

	int rows = 0, cols = 0;
	double* X = nullptr;

	read_matrix_from_disk(file_name, &rows, &cols, &X);

	std::cout << rows << ", " << cols << std::endl;

	free(X);
	X = nullptr;
}

TEST(CWFRTest, hfli) {

	int rows = 0, cols = 0;

	double* X = nullptr;
	double* Y = nullptr;
	double* Z = nullptr;
	double* Sx = nullptr;
	double* Sy = nullptr;

	// oad data
	read_matrix_from_disk("../../data/X.bin", &rows, &cols, &X);
	read_matrix_from_disk("../../data/Y.bin", &rows, &cols, &Y);
	read_matrix_from_disk("../../data/Z.bin", &rows, &cols, &Z);
	read_matrix_from_disk("../../data/Sx.bin", &rows, &cols, &Sx);
	read_matrix_from_disk("../../data/Sy.bin", &rows, &cols, &Sy);

	// map the data to Eigen
	Eigen::Map<MatrixXXd> Xmap(X, rows, cols);
	Eigen::Map<MatrixXXd> Ymap(Y, rows, cols);
	Eigen::Map<MatrixXXd> Zmap(Z, rows, cols);
	Eigen::Map<MatrixXXd> Sxmap(Sx, rows, cols);
	Eigen::Map<MatrixXXd> Symap(Sy, rows, cols);

	std::cout << Zmap.cols() << ", " << Zmap.rows() << std::endl;

	// call the functions
	CWFR wfr(Sxmap, Symap, Xmap, Ymap);
	MatrixXXd Z_calc = wfr();
	MatrixXXd Z_diff = Z_calc - Zmap;

	write_matrix_to_disk("../../data/Z_calc.bin", rows, cols, Z_calc.data());

	free(X);
	free(Y);
	free(Z);
	free(Sx);
	free(Sy);
}