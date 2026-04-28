//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/walks_demo.c: benchmark LAGraph_NumberOfWalks
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Maira Athar, Texas A&M University

//------------------------------------------------------------------------------

// Usage:  WALK_K=4 walks_demo graph.mtx
//         WALK_K=4 walks_demo graph.mtx sources.mtx

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0

// to run with p and p/2 threads, if p = omp_get_max_threads()
// #define NTHREAD_LIST 2
// #define THREAD_LIST 0

// #define NTHREAD_LIST 7
// #define THREAD_LIST 32, 24, 16, 8, 4, 2, 1

// #define NTHREAD_LIST 4
// #define THREAD_LIST 32, 24, 16, 8

// #define NTHREAD_LIST 8
// #define THREAD_LIST 8, 7, 6, 5, 4, 3, 2, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define DEFAULT_K         4
#define N_DEFAULT_SOURCES 64
#define MAX_TRIALS        8

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                     \
{                                                       \
    LAGraph_Delete (&G, msg) ;                          \
    GrB_free (&C) ;                                     \
    GrB_free (&level) ;                                 \
    GrB_free (&parent) ;                                \
    GrB_free (&src_indicator) ;                         \
    GrB_free (&multisrc_vector) ;                       \
    GrB_free (&SourceNodes) ;                           \
    LAGraph_Free ((void **) &t_walks_all,  msg) ;       \
    LAGraph_Free ((void **) &t_msbfs,      msg) ;       \
    LAGraph_Free ((void **) &t_walks_src,  msg) ;       \
    LAGraph_Free ((void **) &t_bfs,        msg) ;       \
}

int main (int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix C = NULL ;
    GrB_Matrix level = NULL ;
    GrB_Vector parent = NULL ;
    GrB_Vector src_indicator  = NULL ;
    GrB_Vector multisrc_vector = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    double *t_walks_all = NULL ;
    double *t_msbfs     = NULL ;
    double *t_walks_src = NULL ;
    double *t_bfs       = NULL ;

    bool burble = false ;
    demo_init (burble) ;

    // read walk length k from WALK_K env var to avoid interfering with
    // readproblem's argv parsing (argv[2] is treated as a sources file)
    int64_t k = DEFAULT_K ;
    const char *k_env = getenv ("WALK_K") ;
    if (k_env != NULL && k_env [0] != '\0') k = atol (k_env) ;
    if (k <= 0) k = DEFAULT_K ;

    //--------------------------------------------------------------------------
    // thread list setup
    //--------------------------------------------------------------------------

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max, nthreads_outer, nthreads_inner ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
    nthreads_max = nthreads_outer * nthreads_inner ;
    if (Nthreads [1] == 0)
    {
        // create thread list automatically
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test:") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_walks_all, nthreads_max+1,
        sizeof (double), msg)) ;
    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_msbfs,     nthreads_max+1,
        sizeof (double), msg)) ;
    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_walks_src, nthreads_max+1,
        sizeof (double), msg)) ;
    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_bfs,       nthreads_max+1,
        sizeof (double), msg)) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, &SourceNodes,
        false, false, false, GrB_INT64, false, argc, argv)) ;

    LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    GrB_Index ntrials ;
    if (SourceNodes == NULL)
    {
        ntrials = (GrB_Index) N_DEFAULT_SOURCES ;
        if (ntrials > n) ntrials = n ;
        GRB_TRY (GrB_Matrix_new (&SourceNodes, GrB_INT64, ntrials, 1)) ;
        for (GrB_Index i = 0 ; i < ntrials ; i++)
        {
            GRB_TRY (GrB_Matrix_setElement_INT64 (SourceNodes,
                (int64_t)(i + 1), i, 0)) ;
        }
    }
    else
    {
        GRB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;
    }
    if (ntrials > MAX_TRIALS) ntrials = MAX_TRIALS ;

    printf ("graph: %s  n: %" PRIu64 "  edges: %" PRIu64 "  k: %" PRId64 "\n",
        matrix_name, (uint64_t) n, (uint64_t) nvals, k) ;
    printf ("source trials: %" PRIu64 "\n", (uint64_t) ntrials) ;
    fflush (stdout) ; fflush (stderr) ;

    //--------------------------------------------------------------------------
    // build MultiSourceBFS source vector
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&multisrc_vector, GrB_INT64, ntrials)) ;
    for (GrB_Index i = 0 ; i < ntrials ; i++)
    {
        int64_t src ;
        GRB_TRY (GrB_Matrix_extractElement_INT64 (&src, SourceNodes, i, 0)) ;
        src-- ;
        GRB_TRY (GrB_Vector_setElement_INT64 (multisrc_vector, src, i)) ;
    }

    //--------------------------------------------------------------------------
    // warmup
    //--------------------------------------------------------------------------

    double twarmup = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, NULL, k)) ;
    twarmup = LAGraph_WallClockTime () - twarmup ;
    GRB_TRY (GrB_free (&C)) ;
    printf ("warmup: NumberOfWalks all-pairs: %g sec\n", twarmup) ;
    fflush (stdout) ; fflush (stderr) ;

    //--------------------------------------------------------------------------
    // benchmark across thread counts
    //--------------------------------------------------------------------------

    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;

        printf ("\n------------------------------- threads: %2d\n", nthreads) ;

        //----------------------------------------------------------------------
        // Part 1: all-pairs NumberOfWalks vs MultiSourceBFS
        //----------------------------------------------------------------------

        {
            double t_run = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, NULL, k)) ;
            t_run = LAGraph_WallClockTime () - t_run ;
            GRB_TRY (GrB_free (&C)) ;
            t_walks_all [nthreads] = t_run ;
            printf ("NumberOfWalks all-pairs  k: %2" PRId64
                "  threads: %2d  time: %10.4f sec\n",
                k, nthreads, t_run) ;
            fflush (stdout) ;
        }

        {
            double t_run = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, NULL, G,
                multisrc_vector, msg)) ;
            t_run = LAGraph_WallClockTime () - t_run ;
            GRB_TRY (GrB_free (&level)) ;
            t_msbfs [nthreads] = t_run ;
            printf ("MultiSourceBFS           k: %2" PRId64
                "  threads: %2d  time: %10.4f sec\n",
                k, nthreads, t_run) ;
            fflush (stdout) ;
        }

        //----------------------------------------------------------------------
        // Part 2: single-source NumberOfWalks vs BFS, averaged over sources
        //----------------------------------------------------------------------

        double total_walks_src = 0, total_bfs = 0 ;
        GRB_TRY (GrB_Vector_new (&src_indicator, GrB_INT64, n)) ;

        for (GrB_Index trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t src ;
            GRB_TRY (GrB_Matrix_extractElement_INT64 (&src, SourceNodes,
                trial, 0)) ;
            src-- ;

            double tb = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGr_BreadthFirstSearch (NULL, &parent, G,
                (GrB_Index) src, msg)) ;
            tb = LAGraph_WallClockTime () - tb ;
            GRB_TRY (GrB_free (&parent)) ;
            total_bfs += tb ;

            GRB_TRY (GrB_Vector_setElement_INT64 (src_indicator, 1,
                (GrB_Index) src)) ;
            double tw = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, src_indicator, k)) ;
            tw = LAGraph_WallClockTime () - tw ;
            GRB_TRY (GrB_free (&C)) ;
            GRB_TRY (GrB_Vector_clear (src_indicator)) ;
            total_walks_src += tw ;

            printf ("trial: %3" PRIu64 "  src: %12" PRId64
                "  walks: %8.4f sec  bfs: %8.4f sec\n",
                (uint64_t) trial, src, tw, tb) ;
            fflush (stdout) ;
        }

        GRB_TRY (GrB_free (&src_indicator)) ;

        t_walks_src [nthreads] = total_walks_src / (double) ntrials ;
        t_bfs       [nthreads] = total_bfs       / (double) ntrials ;

        //----------------------------------------------------------------------
        // summary (printed to stderr for script parsing)
        //----------------------------------------------------------------------

        printf (         "Avg: NumberOfWalks all-pairs    k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_all [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: NumberOfWalks all-pairs    k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_all [nthreads], matrix_name) ;

        printf (         "Avg: MultiSourceBFS             k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_msbfs [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: MultiSourceBFS             k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_msbfs [nthreads], matrix_name) ;

        printf (         "Avg: NumberOfWalks single-src   k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_src [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: NumberOfWalks single-src   k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_src [nthreads], matrix_name) ;

        printf (         "Avg: BFS single-source          k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_bfs [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: BFS single-source          k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_bfs [nthreads], matrix_name) ;

        fflush (stdout) ; fflush (stderr) ;
    }

    // restore default thread count
    LAGRAPH_TRY (LAGraph_SetNumThreads (nthreads_outer, nthreads_inner, msg)) ;

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
