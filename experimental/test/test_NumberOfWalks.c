//----------------------------------------------------------------------------
// LAGraph/experimental/test/test_NumberOfWalks.c: test for LAGraph_NumberOfWalks
//----------------------------------------------------------------------------

// Test for LAGraph_NumberOfWalks: Recursive binary exponentiation for 
// computing A^k (number of walks of length k)

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

char msg [LAGRAPH_MSG_LEN] ;

void test_NumberOfWalks_path_graph (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, C = NULL ;

    // test with the path graph 0-1-2

    // create adjacency matrix for path: 0-1-2
    GrB_Index n = 3;
    OK (GrB_Matrix_new (&A, GrB_INT64, n, n)) ;

    GrB_Index i[4] = {0, 1, 1, 2};
    GrB_Index j[4] = {1, 0, 2, 1};
    int64_t x[4] = {1, 1, 1, 1};
    OK (GrB_Matrix_build (A, i, j, x, 4, GrB_PLUS_INT64)) ;

    // compute A^2 (walks of length 2)
    OK (LAGraph_NumberOfWalks (&C, A, 2)) ;

    // check result: 
    // C(1,1) -> 2
    // C(0,2) -> 1
    int64_t val11, val02;
    OK (GrB_Matrix_extractElement_INT64 (&val11, C, 1, 1)) ;
    OK (GrB_Matrix_extractElement_INT64 (&val02, C, 0, 2)) ;

    printf ("\nA^2 for path graph 0-1-2:\n") ;
    printf ("C(1,1) = %lld (expected 2)\n", (long long) val11) ;
    printf ("C(0,2) = %lld (expected 1)\n", (long long) val02) ;

    TEST_CHECK (val11 == 2) ;
    TEST_CHECK (val02 == 1) ;

    // free everything and finalize LAGraph
    OK (GrB_free (&A)) ;
    OK (GrB_free (&C)) ;

    LAGraph_Finalize (msg) ;
}


TEST_LIST =
{
    {"NumberOfWalks_path_graph", test_NumberOfWalks_path_graph},
    {NULL, NULL}
} ;
