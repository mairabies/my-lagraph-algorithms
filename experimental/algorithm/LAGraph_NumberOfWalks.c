#include "GraphBLAS.h"
#include "LAGraphX.h"

/**
 * NumberOfWalks: Recursive implementation using binary exponentiation.
 * 
 * @param C      The output matrix (A^k)
 * @param A      The adjacency matrix
 * @param k      Walk length
 */

GrB_Info LAGraph_NumberOfWalks (GrB_Matrix *C, GrB_Matrix A, int64_t k) 
{
    if (C == NULL || A == NULL || k < 0) return GrB_INVALID_VALUE;

    GrB_Info info;
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);

    // BASE CASES: k=0, 1, 2

    // identity matrix, is a for loop necessary here?
    if (k == 0) {
        GrB_Matrix_new(C, GrB_INT64, n, n);
        for (GrB_Index i = 0; i < n; i++) GrB_Matrix_setElement_INT64(*C, 1, i, i);
        return GrB_SUCCESS;
    }

    // adjacency matrix itself
    if (k == 1) return GrB_Matrix_dup(C, A);
    
    // A^2
    if (k == 2) {
        GrB_Matrix_new(C, GrB_INT64, n, n);
        return GrB_mxm(*C, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, A, A, NULL);
    }

    // RECURSION
    GrB_Matrix T;
    info = LAGr_NumberOfWalks(&T, A, k / 2);
    if (info != GrB_SUCCESS) return info;

    // square res = T^2
    // Res is to hold the result of T^2, is this necessary?
    GrB_Matrix Res;
    GrB_Matrix_new(&Res, GrB_INT64, n, n);
    GrB_mxm(Res, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, T, T, NULL);
    GrB_free(&T);

    // If k is odd, multiply by A (Res = Res * A)
    if (k % 2 != 0) {
        GrB_Matrix_new(C, GrB_INT64, n, n);
        GrB_mxm(*C, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, Res, A, NULL);
        GrB_free(&Res);
    // if k is even, just return T^2 (Res)
    } else {
        *C = Res;
    }

    return GrB_SUCCESS;
}