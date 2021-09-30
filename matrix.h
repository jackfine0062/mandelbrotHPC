/**
 * Declares all of the matrix functions (which are defined in matrix.c).
 */

#pragma once

#ifndef __STDC_WANT_LIB_EXT2__
#define __STDC_WANT_LIB_EXT2__ 1 // allows some extra features in C
#endif
#ifndef _XOPEN_SOURCE
#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif
#endif
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

struct _Matrix {
    // Our basic matrix structure
    size_t rows, cols, size; // size is simply rows*cols, but it comes up a lot
    float* data;
    char data_source; // one of DATA_MALLOCED, DATA_MEMMAPPED, or DATA_BORROWED
};
typedef struct _Matrix Matrix; // make type "struct _Matrix" just "Matrix"

typedef float (*unary_func)(float);
typedef float (*binary_func)(float, float);


#ifdef __cplusplus
extern "C" {
#endif

//////////////////// Matrix Creation Functions //////////////////// 

/**
 * Creates a new matrix of the given rows and columns. The data is NOT
 * initialized. Returns the newly created matrix. The data_source attribute is
 * set to DATA_MALLOCED.
 */
Matrix* matrix_create_raw(size_t rows, size_t cols);

/**
 * Frees a Matrix object. Depending on the data_source, either the data is
 * free()ed, munmap()ed, or nothing is done to it. Afterwards the Matrix
 * variable itself is freed.
 */
void matrix_free(Matrix* M);

/**
 * Creates a new matrix of the given rows and columns. The data is all set to
 * zeroes. Returns the newly created matrix.
 */
Matrix* matrix_zeros(size_t rows, size_t cols);

/**
 * Fills a matrix with all 0s.
 */
void matrix_fill_zeros(Matrix* M);

/**
 * Creates a new matrix of the given rows and columns. The data is all set to
 * zeroes except the main diagonal which is set to ones. Returns the newly
 * created matrix.
 */
Matrix* matrix_identity(size_t rows, size_t cols);

/**
 * Creates a new matrix of the given rows and columns. The data is filled in
 * random values in [0.0, 1.0]. Returns the newly created matrix.
 * 
 * The caller is responsible for making sure srand() is appropriately called
 * before this function is called.
 */
Matrix* matrix_random(size_t rows, size_t cols);

/**
 * Creates a new matrix which is a copy of the given matrix.
 */
Matrix* matrix_copy(const Matrix* M);

/**
 * Creates a new matrix by loading the data from the given CSV file. It is
 * assumed that every row in the file has the same number of values. If any
 * row has more values than the first row, those extra values are ignored. If
 * any row has less than the first row, the missing data is filled with 0s. If
 * there is a problem reading from the file, NULL is returned.
 */
Matrix* matrix_from_csv(FILE* file);

/**
 * Same as matrix_from_csv() but takes a file path instead.
 */
Matrix* matrix_from_csv_path(const char* path);

/**
 * Saves a matrix to a CSV file.
 * 
 * If the file argument is given as stdout, this will print it to the terminal.
 */
void matrix_to_csv(FILE* file, const Matrix* M);

/**
 * Same as matrix_to_csv() but takes a file path instead.
 */
bool matrix_to_csv_path(const char* path, const Matrix* M);

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
Matrix* matrix_from_npy(FILE* file);

/**
 * Same as matrix_from_npy() but takes a file path instead.
 */
Matrix* matrix_from_npy_path(const char* path);

/**
 * Saves a matrix to a NPY file. This is a file format used by the numpy
 * library. This will return false if the data cannot be written.
 */
bool matrix_to_npy(FILE* file, const Matrix* M);

/**
 * Same as matrix_to_npy() but takes a file path instead.
 */
bool matrix_to_npy_path(const char* path, const Matrix* M);


//////////////////// Matrix Comparison Functions //////////////////// 

/**
 * Check if the two matrices have the same rows and columns and that all values
 * in two matrices are exactly equal. Note that due to floating-point math
 * inaccuracies, this will frequently be false even when performing the same
 * exact math. Use matrix_allclose() instead.
 */
bool matrix_equal(const Matrix* A, const Matrix* B);

/**
 * Check if all values in two matrices are within a relative and absolute
 * tolerance of each other. This returns true if, for each pair of elements in
 * A and B, abs(a-b) <= (atol+rtol*abs(b)) is true. Good values for rtol and
 * atol are 1e-05 and 1e-08. If any value in either matrix is nan, this will
 * always compare as false.
 */
bool matrix_allclose(const Matrix* A, const Matrix* B, float rtol, float atol);


//////////////////// Basic Matrix Functions //////////////////// 
// NOTES:
//  * all of the *_apply() functions operate in-place, modifying the first
//    matrix argument
//  * all of the *_map() functions create a new output matrix and return it
//  * everything before the apply or map tells you the argument types (either
//    a Matrix* or a float [for scalar])
//  * all of these functions are incredibly similar, just a few changes between
//    each one


/**
 * Apply a unary function to every element of a matrix replacing the value with
 * the return value of the function. This operates in-place.
 * 
 * For example, the unary function could be one of the built-in functions like
 * sin, cos, fabs, exp, log, sqrt, floor, or ceil. It will also work with any
 * custom functions that have a single float parameter and return a float.
 */
void matrix_apply(Matrix* M, unary_func func);

/**
 * Apply a unary function to every element of a matrix and save the value to a
 * new matrix.
 * 
 * For example, the unary function could be one of the built-in functions like
 * sin, cos, fabs, exp, log, sqrt, floor, or ceil. It will also work with any
 * custom functions that have a single float parameter and return a float.
 */
Matrix* matrix_map(const Matrix* M, unary_func func);

/**
 * Apply a binary function to every pair of elements from two matrices,
 * replacing the value in the first matrix with the return value of the
 * function. This operates in-place.
 * 
 * This returns false if the two matrices are not the same shape.
 * 
 * For example, the binary function could be one of the built-in functions like
 * pow, fmod, fmax, fmin, hypot. It will also work with any custom functions that
 * have two float paramters and return a float value.
 */
bool matrix_matrix_apply(Matrix* A, const Matrix* B, binary_func func);

/**
 * Apply a binary function to every pair of elements from two matrices, and
 * save the return value to a new matrix.
 * 
 * This returns NULL if the two matrices are not the same shape.
 * 
 * For example, the binary function could be one of the built-in functions like
 * pow, fmod, fmax, fmin, hypot. It will also work with any custom functions that
 * have two float paramters and return a float value.
 */
Matrix* matrix_matrix_map(const Matrix* A, const Matrix* B, binary_func func);

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the second argument, replacing the value in the matrix
 * with the return value of the function. This operates in-place.
 * 
 * For example, the binary function could be one of the built-in functions like
 * pow, fmod, fmax, fmin, hypot. It will also work with any custom functions that
 * have two float paramters and return a float value.
 */
void matrix_scalar_apply(Matrix* A, float b, binary_func func);

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the second argument, saving the return value to a new
 * matrix.
 * 
 * For example, the binary function could be one of the built-in functions like
 * pow, fmod, fmax, fmin, hypot. It will also work with any custom functions that
 * have two float paramters and return a float value.
 */
Matrix* matrix_scalar_map(const Matrix* A, float b, binary_func func);

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the first argument, replacing the value in the matrix
 * with the return value of the function. This operates in-place.
 * 
 * For example, the binary function could be one of the built-in functions like
 * pow, fmod, fmax, fmin, hypot. It will also work with any custom functions that
 * have two float paramters and return a float value.
 */
void scalar_matrix_apply(float a, Matrix* B, binary_func func);

/**
 * Apply a binary function to every element in the matrix with the scalar
 * (i.e. fixed value) as the first argument, saving the return value to a new
 * matrix.
 * 
 * For example, the binary function could be one of the built-in functions like
 * pow, fmod, fmax, fmin, hypot. It will also work with any custom functions that
 * have two float paramters and return a float value.
 */
Matrix* scalar_matrix_map(float a, const Matrix* B, binary_func func);

/**
 * Apply a binary function to every elements and the previous return value,
 * reducing down to a single value. Useful for things like fmax().
 */
float matrix_reduce(const Matrix* A, binary_func func);

/**
 * Apply a binary function to every elements and the previous return value,
 * reducing down to a single value for each row. Useful for things like fmax().
 */
Matrix* matrix_reduce_rows(const Matrix* A, binary_func func);

/**
 * Apply a binary function to every elements and the previous return value,
 * reducing down to a single value for each column. Useful for things like
 * fmax().
 */
Matrix* matrix_reduce_cols(const Matrix* A, binary_func func);


//////////////////// Unary and Binary Functions //////////////////// 
// These are to be used with the *_apply() and *_map() functions like:
//   M = matrix_map(A, sin)    // creates a new matrix where each value is the sin of the corresponding value in A
//   matrix_apply(M, negative) // replaces all values in M with their negative
//   C = matrix_matrix_map(A, B, subtract) // all values in C are the difference of values in A and B

// Unary functions
static inline float reciprocal(float a) { return 1 / a; }
static inline float positive(float a) { return +a; }
static inline float negative(float a) { return -a; }
// included in math.h:
// fabs, exp, exp2, expm1, log, log10, log2, log1p, sqrt, cbrt
// sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh
// erf, erfc, tgamma, lgamma
// ceil, floor, trunc, round, nearbyint, rint
// and more

// Binary functions
static inline float add(float a, float b) { return a + b; }
static inline float subtract(float a, float b) { return a - b; }
static inline float multiply(float a, float b) { return a * b; }
static inline float divide(float a, float b) { return a / b; }
// included in math.h:
// fmod, remainder, fmin, fmax, fdim, pow, hypot, atan2


//////////////////// Matrix Functions ////////////////////

/**
 * MATRIX_AT() works like matrix_get() and matrix_set() combined from homework 1
 * except that it doesn't require an actual function call and can be used on
 * either the right or left side of an = (see matrix_multiplication for an
 * example).
 */
#define MATRIX_AT(M, i, j) ((M)->data[(i)*(M)->cols+(j)])

/**
 * Matrix multiplication. https://en.wikipedia.org/wiki/Matrix_multiplication
 * Output is written to the third argument which must be the right size.
 * If the inner dimensions are not equal or output not the right size, this
 * returns false.
 */
bool matrix_multiplication(const Matrix* A, const Matrix* B, Matrix* C);

#ifdef __cplusplus
}
#endif
