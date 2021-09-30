/**
 * Matrix definitions
 */

#include "matrix.h"
#include "matrix_io_helpers.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>

#include <unistd.h>
#include <sys/mman.h>

#define DATA_MALLOCED   1 // data is from malloc(), needs free()
#define DATA_MEMMAPPED  2 // data is from mmap(), needs munmap()
#define DATA_BORROWED   3 // data is from elsewhere and should not be freed


//////////////////// Matrix Creation Functions //////////////////// 

/**
 * Creates a Matrix with attributes set with the arguments.
 */
Matrix* matrix_alloc(size_t rows, size_t cols, float* data, char data_source) {
    Matrix* M = (Matrix*)malloc(sizeof(Matrix));
    M->rows = rows;
    M->cols = cols;
    M->size = rows * cols;
    M->data = data;
    M->data_source = data_source;
    return M;
}

/**
 * Creates a new matrix of the given rows and columns. The data is NOT
 * initialized. Returns the newly created matrix.
 */
Matrix* matrix_create_raw(size_t rows, size_t cols) {
    return matrix_alloc(rows, cols,
        (float*)malloc(rows*cols*sizeof(float)), DATA_MALLOCED);
}

/**
 * Frees a Matrix object and its data (if not borrowed).
 */
void matrix_free(Matrix* M) {
    if (M->data_source == DATA_MEMMAPPED) {
        size_t addr = ((size_t)M->data) & ~(sysconf(_SC_PAGE_SIZE)-1);
        munmap((void*)addr, M->size*sizeof(float));
    } else if (M->data_source == DATA_MALLOCED) {
        free(M->data);
    }
    free(M);
}

/**
 * Creates a new matrix of the given rows and columns. The data is all set to
 * zeroes. Returns the newly created matrix.
 */
Matrix* matrix_zeros(size_t rows, size_t cols) {
    Matrix* M = matrix_create_raw(rows, cols);
    memset(M->data, 0, M->size*sizeof(float));
    return M;
}

/**
 * Fills a matrix with all 0s.
 */
void matrix_fill_zeros(Matrix* M) {
    memset(M->data, 0, M->size*sizeof(float));
}

/**
 * Creates a new matrix of the given rows and columns. The data is all set to
 * zeroes except the main diagonal which is set to ones. Returns the newly
 * created matrix.
 */
Matrix* matrix_identity(size_t rows, size_t cols) {
    Matrix* M = matrix_zeros(rows, cols);
    size_t size = rows < cols ? rows : cols;
    for (size_t i = 0, index = 0; i < size; i++, index += cols+1) {
        M->data[index] = 1.0;
    }
    return M;
}

/**
 * Creates a new matrix of the given rows and columns. The data is filled in
 * random values in [0.0, 1.0]. Returns the newly created matrix.
 * 
 * The caller is responsible for making sure srand() is appropriately called
 * before this function is called.
 */
Matrix* matrix_random(size_t rows, size_t cols) {
    Matrix* M = matrix_create_raw(rows, cols);
    for (size_t i = 0; i < M->size; i++) {
        M->data[i] = (float)rand() / RAND_MAX;
    }
    return M;
}

/**
 * Creates a new matrix which is a copy of the given matrix.
 */
Matrix* matrix_copy(const Matrix* M) {
    Matrix* out = matrix_create_raw(M->rows, M->cols);
    memcpy(out->data, M->data, M->size*sizeof(float));
    return out;
}

/**
 * Creates a new matrix by loading the data from the given CSV file. It is
 * assumed that every row in the file has the same number of values. If any
 * row has more values than the first row, those extra values are ignored. If
 * any row has less than the first row, the missing data is filled with 0s. If
 * there is a problem reading from the file, NULL is returned.
 */
Matrix* matrix_from_csv(FILE* file) {
    // Get the first line from the file
    char* line = NULL;
    size_t len = 0;
    if (getline(&line, &len, file) <= 0 || len <= 0) {
        free(line);
        return NULL;
    }

    // Parse the first line of the file, getting the number of columns
    float* data = NULL;
    size_t rows = 1, cols = __read_csv_first_line(line, &data);
    if (cols == 0) { free(line); free(data); return NULL; }

    // Get remaining rows of the file
    size_t rows_alloc = cols; // assume square, adjust later if needed
    data = (float*)realloc(data, cols*rows_alloc*sizeof(float));
    while (!feof(file)) {
        if (getline(&line, &len, file) <= 0) { break; }
        __read_csv_line(line, &(data[rows*cols]), cols);
        if (++rows == rows_alloc) {
            rows_alloc += rows_alloc > 128 ? 128 : rows_alloc;
            data = (float*)realloc(data, rows_alloc*cols*sizeof(float));
        }
    }
    free(line);

    // Make the final matrix
    return matrix_alloc(rows, cols,
        (float*)realloc(data, rows*cols*sizeof(float)), DATA_MALLOCED);
}

/**
 * Same as matrix_from_csv() but takes a file path instead.
 */
Matrix* matrix_from_csv_path(const char* path) {
    FILE* f = fopen(path, "r+");
    if (!f) { return NULL; }
    Matrix* M = matrix_from_csv(f);
    fclose(f);
    return M;
}

/**
 * Saves a matrix to a CSV file.
 * 
 * If the file argument is given as stdout, this will print it to the terminal.
 */
void matrix_to_csv(FILE* file, const Matrix* M) {
    if (M->rows < 1 || M->cols < 1) { return; }
    for (size_t i = 0, ind = 0; i < M->rows; i++, ind++) {
        for (size_t j = 0; j < M->cols-1; j++, ind++) {
            fprintf(file, "%f,", M->data[ind]);
        }
        fprintf(file, "%f\n", M->data[ind]);
    }
}

/**
 * Same as matrix_to_csv() but takes a file path instead.
 */
bool matrix_to_csv_path(const char* path, const Matrix* M) {
    FILE* f = fopen(path, "w");
    if (!f) { return false; }
    matrix_to_csv(f, M);
    fclose(f);
    return true;
}

/**
 * Creates a new matrix by loading the data from the given NPY file. This is
 * a file format used by the numpy library. This function only supports arrays
 * that are little-endian floats, c-contiguous, and 1 or 2 dimensional. The
 * file is loaded as memory-mapped so it is backed by the file and loaded
 * on-demand. The file should be opened for reading or reading and writing.
 * 
 * This will return NULL if the data cannot be read, the file format is not
 * recognized, there are memory allocation issues, or the array is not a
 * supported shape or data type.
 */
Matrix* matrix_from_npy(FILE* file) {
    // Read the header, check it, and get the shape of the matrix
    size_t sh[2], offset;
    if (!__npy_read_header(file, sh, &offset)) { return NULL; }

    // Get the memory mapped data
    void* x = (void*)mmap(NULL, sh[0]*sh[1]*sizeof(float) + offset,
                          PROT_READ|PROT_WRITE, MAP_SHARED, fileno(file), 0);
    if (x == MAP_FAILED) { return NULL; }

    // Make the matrix itself
    float* data = (float*)(((char*)x) + offset);
    return matrix_alloc(sh[0], sh[1], data, DATA_MEMMAPPED);
}

/**
 * Same as matrix_from_npy() but takes a file path instead.
 */
Matrix* matrix_from_npy_path(const char* path) {
    FILE* f = fopen(path, "r+b");
    if (!f) { return NULL; }
    Matrix* M = matrix_from_npy(f);
    fclose(f);
    return M;
}

/**
 * Saves a matrix to a NPY file. This is a file format used by the numpy
 * library. This will return false if the data cannot be written.
 */
bool matrix_to_npy(FILE* file, const Matrix* M) {
    // create the header
    char header[128];
    size_t len = snprintf(header, sizeof(header), "\x93NUMPY\x01   "
        "{'descr': '<f4', 'fortran_order': False, 'shape': (%zu, %zu), }",
        M->rows, M->cols);
    if (len < 0) { return false; }
    header[7] = 0; // have to after the string is written
    *(unsigned short*)&header[8] = sizeof(header) - 10;
    memset(header + len, ' ', sizeof(header)-len-1);
    header[sizeof(header)-1] = '\n';

    // write the header and the data
    return fwrite(header, 1, sizeof(header), file) == sizeof(header) &&
        fwrite(M->data, sizeof(float), M->size, file) == M->size;
}

/**
 * Same as matrix_to_npy() but takes a file path instead.
 */
bool matrix_to_npy_path(const char* path, const Matrix* M) {
    FILE* f = fopen(path, "wb");
    if (!f) { return false; }
    matrix_to_npy(f, M);
    fclose(f);
    return true;
}


//////////////////// Matrix Comparison Functions //////////////////// 

/**
 * Check if the two matrices have the same rows and columns and that all values
 * in two matrices are exactly equal. Note that due to floating-point math
 * inaccuracies, this will frequently be false even when performing the same
 * exact math. Use matrix_allclose() instead.
 */
bool matrix_equal(const Matrix* A, const Matrix* B) {
    return A->rows == B->rows && A->cols == B->cols &&
        memcmp(A->data, B->data, A->size*sizeof(float)) == 0;
}

/**
 * Check if all values in two matrices are within a relative and absolute
 * tolerance of each other. This returns true if, for each pair of elements in
 * A and B, abs(a-b) <= (atol+rtol*abs(b)) is true. Good values for rtol and
 * atol are 1e-05 and 1e-08. If any value in either matrix is nan, this will
 * always compare as false.
 */
bool matrix_allclose(const Matrix* A, const Matrix* B, float rtol, float atol) {
    if (A->rows != B->rows || A->cols != B->cols) { return false; }
    for (size_t i = 0; i < A->size; i++) {
        if (fabs(A->data[i] - B->data[i]) >= atol + rtol * fabs(B->data[i])) {
            return false;
        }
    }
    return true;
}


//////////////////// Basic Matrix Functions //////////////////// 

/**
 * Apply a unary function to every element of a matrix replacing the value with
 * the return value of the function. This operates in-place.
 */
void matrix_apply(Matrix* M, unary_func func) {
    for (size_t i = 0; i < M->size; i++) { M->data[i] = func(M->data[i]); }
}

/**
 * Apply a unary function to every element of a matrix and save the value to a
 * new matrix.
 */
Matrix* matrix_map(const Matrix* M, unary_func func) {
    Matrix* out = matrix_create_raw(M->rows, M->cols);
    for (size_t i = 0; i < M->size; i++) { out->data[i] = func(M->data[i]); }
    return out;
}

/**
 * Apply a binary function to every pair of elements from two matrices,
 * replacing the value in the first matrix with the return value of the
 * function. This operates in-place. This returns false if the two matrices are
 * not the same shape.
 */
bool matrix_matrix_apply(Matrix* A, const Matrix* B, binary_func func) {
    if (A->rows != B->rows || A->cols != B->cols) { return false; }
    for (size_t i = 0; i < A->size; i++) {
        A->data[i] = func(A->data[i], B->data[i]);
    }
    return true;
}

/**
 * Apply a binary function to every pair of elements from two matrices, and
 * save the return value to a new matrix. This returns NULL if the two matrices
 * are not the same shape.
 */
Matrix* matrix_matrix_map(const Matrix* A, const Matrix* B, binary_func func) {
    if (A->rows != B->rows || A->cols != B->cols) { return NULL; }
    Matrix* out = matrix_create_raw(A->rows, A->cols);
    for (size_t i = 0; i < A->size; i++) {
        out->data[i] = func(A->data[i], B->data[i]);
    }
    return out;
}

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the second argument, replacing the value in the matrix
 * with the return value of the function. This operates in-place.
 */
void matrix_scalar_apply(Matrix* A, float b, binary_func func) {
    for (size_t i = 0; i < A->size; i++) { A->data[i] = func(A->data[i], b); }
}

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the second argument, saving the return value to a new
 * matrix.
 */
Matrix* matrix_scalar_map(const Matrix* A, float b, binary_func func) {
    Matrix* out = matrix_create_raw(A->rows, A->cols);
    for (size_t i = 0; i < A->size; i++) { out->data[i] = func(A->data[i], b); }
    return out;
}

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the first argument, replacing the value in the matrix
 * with the return value of the function. This operates in-place.
 */
void scalar_matrix_apply(float a, Matrix* B, binary_func func) {
    for (size_t i = 0; i < B->size; i++) { B->data[i] = func(a, B->data[i]); }
}

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the first argument, saving the return value to a new
 * matrix.
 */
Matrix* scalar_matrix_map(float a, const Matrix* B, binary_func func) {
    Matrix* out = matrix_create_raw(B->rows, B->cols);
    for (size_t i = 0; i < B->size; i++) { out->data[i] = func(a, B->data[i]); }
    return out;
}

/**
 * Apply a binary function to every elements and the previous return value,
 * reducing down to a single value. Useful for things like fmax().
 */
float matrix_reduce(const Matrix* A, binary_func func) {
    float val = A->data[0];
    for (size_t i = 1; i < A->size; i++) { val = func(A->data[i], val); }
    return val;
}

/**
 * Apply a binary function to every elements and the previous return value,
 * reducing down to a single value for each row. Useful for things like fmax().
 */
Matrix* matrix_reduce_rows(const Matrix* A, binary_func func) {
    Matrix* out = matrix_create_raw(A->rows, 1);
    for (size_t i = 0; i < A->rows; i++) {
        float val = MATRIX_AT(A, i, 0);
        for (size_t j = 1; j < A->cols; j++) {
            val = func(MATRIX_AT(A, i, j), val);
        }
        MATRIX_AT(out, i, 0) = val;
    }
    return out;
}

/**
 * Apply a binary function to every elements and the previous return value,
 * reducing down to a single value for each column. Useful for things like
 * fmax().
 */
Matrix* matrix_reduce_cols(const Matrix* A, binary_func func) {
    Matrix* out = matrix_create_raw(1, A->cols);
    for (size_t i = 0; i < A->cols; i++) {
        float val = MATRIX_AT(A, 0, i);
        for (size_t j = 1; j < A->rows; j++) {
            val = func(MATRIX_AT(A, j, i), val);
        }
        MATRIX_AT(out, 0, i) = val;
    }
    return out;
}


//////////////////// Matrix Functions ////////////////////

/**
 * Matrix multiplication. https://en.wikipedia.org/wiki/Matrix_multiplication
 * Output is written to the third argument which must be the right size.
 * If the inner dimensions are not equal or output not the right size, this
 * returns false.
 */
bool matrix_multiplication(const Matrix* A, const Matrix* B, Matrix* C) {
    const size_t m = A->rows, n = A->cols, p = B->cols;
    if (n != B->rows || m != C->rows || p != C->cols) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < p; j++) {
            float value = 0;
            for (size_t k = 0; k < n; k++) {
                value += MATRIX_AT(A, i, k) * MATRIX_AT(B, k, j);
            }
            MATRIX_AT(C, i, j) = value;
        }
    }
    
    return true;
}
