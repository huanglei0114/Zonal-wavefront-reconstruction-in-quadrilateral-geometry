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
	print_matrix_in_matlab_format(rows, cols, X);
}