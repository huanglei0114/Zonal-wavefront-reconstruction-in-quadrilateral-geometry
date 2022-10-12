#ifndef MATRIX_IO
#define MATRIX_IO

#include <cstdio>
#include <cstdlib>
#include <cctype>

inline int ID_1D(int x, int y, int width) { return (y * width + x); }

template <class T>
int read_matrix_from_disk(const char* filename, int* rows, int* cols, T** matrix)
{
	FILE* file;
	fopen_s(&file, filename, "rb");
	if (!file)
	{
		printf("Can't open input matrix file: %s.\n", filename);
		return 1;
	}

	if (read_matrix_size_from_stream(file, rows, cols) != 0)
	{
		printf("Error reading matrix header from disk file: %s.\n", filename);
		return 1;
	}

	//int size = (*m) * (*n) * sizeof(T) + 2 * sizeof(int);
	*matrix = (T*)malloc(sizeof(T) * (*rows) * (*cols));

	if (read_matrix_from_stream(file, *rows, *cols, *matrix) != 0)
	{
		printf("Error reading matrix data from disk file: %s.\n", filename);
		return 1;
	}

	fclose(file);

	return 0;
}

template <class T>
int read_matrix_from_stream(FILE* file, int rows, int cols, T* matrix)
{
	unsigned int readBytes;
	if ((readBytes = fread(matrix, sizeof(T), rows * cols, file)) < (unsigned int)rows * cols)
	{
		printf("Error: I have only read %u bytes. sizeof(T)=%lu\n", readBytes, sizeof(T));
		return 1;
	}

	return 0;
}

template <class T>
int write_matrix_to_disk(const char* filename, int rows, int cols, T* matrix)
{
	FILE* file;
	fopen_s(&file, filename, "wb");
	if (!file)
	{
		printf("Can't open output file: %s.\n", filename);
		return 1;
	}

	if (write_matrix_header_to_stream(file, rows, cols) != 0)
	{
		printf("Error writing the matrix header to disk file: %s.\n", filename);
		return 1;
	}

	if (write_matrix_to_stream(file, rows, cols, matrix) != 0)
	{
		printf("Error writing the matrix to disk file: %s.\n", filename);
		return 1;
	}

	fclose(file);

	return 0;
}

template <class T>
int write_matrix_to_stream(FILE* file, int rows, int cols, T* matrix)
{
	if ((int)(fwrite(matrix, sizeof(T), rows * cols, file)) < rows * cols)
		return 1;
	return 0;
}

template <class T>
void print_matrix_in_matlab_format(int rows, int cols, T* U)
{
	printf("{\n");

	for (int i = 0; i < rows; i++)
	{
		printf("{");

		for (int j = 0; j < cols; j++)
		{
			printf("%f", U[ID_1D(j, i, cols)]);
			if (j != cols - 1)
				printf(", ");
		}
		printf("}");

		if (i != rows - 1)
			printf(",\n");
	}

	printf("\n}\n");
}

static int read_matrix_size_from_stream(FILE* file, int* rows, int* cols)
{
	if (fread(rows, sizeof(int), 1, file) < 1)
		return 1;
	if (fread(cols, sizeof(int), 1, file) < 1)
		return 1;

	return 0;
}

static int write_matrix_header_to_stream(FILE* file, int rows, int cols)
{
	if (fwrite(&rows, sizeof(int), 1, file) < 1)
		return 1;
	if (fwrite(&cols, sizeof(int), 1, file) < 1)
		return 1;

	return 0;
}


#endif // !MATRIX_IO
