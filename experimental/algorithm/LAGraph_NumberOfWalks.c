#include <stdio.h>
#include <stdlib.h>
#include "GraphBLAS.h"
#include "LAGraphX.h"

/**
 * simple brute force imlplemetation for now
 * O(k) matrix multiplications where k is walk length from i to j
 * 
 * May optimize with elimiating repeated squaring. 
 * Looking into binary exponentiation
**/

GrB_Info NumberOfWalks (GrB_Matrix *C, GrB_Matrix A, int64_t k) 
{
    GrB_Info info;
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);

    // case k=0 - the number of walks of length 0 is the identity matrix
    if (k == 0) {
        GrB_Matrix_new(C, GrB_INT64, n, n);
        for (GrB_Index i = 0; i < n; i++) {
            GrB_Matrix_setElement_INT64(*C, 1, i, i);
        }
        return GrB_SUCCESS;
    }

    // initialize result as a copy of A
    GrB_Matrix_dup(C, A);

    // walk length 1 is just the adjacency matrix itself, so we can return early
    if (k == 1) return GrB_SUCCESS;  

    //  multiply by A, k-1 more times
    for (int64_t i = 1; i < k; i++) 
    {
        GrB_Matrix Temp;
        GrB_Matrix_new(&Temp, GrB_INT64, n, n);

        // temp = C * A
        // PLUS_TIMES semiring (sum of products)
        GrB_mxm(Temp, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, *C, A, NULL);

        // free old result and update it with the new one
        GrB_free(C);
        *C = Temp;
    }
    return GrB_SUCCESS;
}