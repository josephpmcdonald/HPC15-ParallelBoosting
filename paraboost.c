#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"
#include "util.h"
#include "mpi.h"


/* paraboost.c implements a parallelized version of AdaBoost, kept in a
 * separate file to simplify compilation.
 */


double WeakLearner(Node *tree, double *x) {

/* WeakLearner is a function that takes an observation/vector of features and
 * returns a prediction on the label of that observation. The parameter tree is
 * a pointer to the root of the decision tree that WeakLearner represents as a
 * base classifier.
 *
 * Duplicate of TestPoint.
 *
 */

    int ind;
    Node *node = tree;

    while (node->left != NULL) {
        ind = node->index;
        if (x[ind] <= node->threshold)
            node = node->left;
        else
            node = node->right;
    }

    return node->label;
}


double PError(Node *tree, double **data, Pod **base, int n) {

/* Note that Error returns the weighted error on the data.
 *
 */

    int i;
    double error = 0;

    for (i = 0; i < n; ++i) {
        if (WeakLearner(tree, data[i])*base[i]->label < 0)
            error += base[i]->weight;
    }

    return error;
}


int main (int argc, char *argv[]) {

/*
 */

    MPI_Init(&argc, &argv);
    if (argc != 2) {
        printf("Must supply T as second argument\n");
        abort();
    }
    int T = atoi(argv[1]);
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (p > D-1) {
        printf("Too many processors. Aborting...\n");
        abort();
    }

    ///////SCATTERING///THE///DATA//////////////////////////////////
    int num_features = (D-1)/p;
    int remainder = (D-1)%p;
    int n = N;
    int i, j;
    int *feature_list;
    if (rank < remainder) {
        num_features += 1;
        feature_list = malloc(num_features*sizeof(int));
        for (i = 0; i < num_features; ++i)
            feature_list[i] = rank*num_features + i;
    }
    else {
        feature_list = malloc(num_features*sizeof(int));
        for (i = 0; i < num_features; ++i)
            feature_list[i] = remainder*(num_features+1) + (rank - remainder)*num_features + i;
    }
    printf("Rank %d has feature %d through %d\n", rank, feature_list[0], feature_list[num_features-1]);
    

    //////COPYDATA///////////////////////////////////////////
    double **MYDATA = ParHIGGS(feature_list, num_features);

    return;

    //Allocate pointers to pods 
    Pod **base = malloc(n*sizeof(Pod*));

    //Each process stores array of pods containing feature
    for (i = 0; i < n; ++i) {
        base[i] = malloc(sizeof(Pod));
        base[i]->val = malloc(num_features*sizeof(double));
    }

    for (i = 0; i < n; ++i) {
        base[i]->key = i; //key
        for (j = 0; j < num_features; j++)
            base[i]->val[j] = MYDATA[i][j]; //feature
        base[i]->label = MYDATA[i][num_features]; //label
        base[i]->weight = 1./n; //weight
    }

    Pod ***data = malloc(num_features*sizeof(Pod**));
    int feat;
    for (feat = 0; feat < num_features; feat++){
        data[feat] = malloc(n*sizeof(Pod*));
        for (i = 0; i < n; ++i)
            data[feat][i] = base[i];
    }

    free(MYDATA);

    double **ALLDATA = NULL;
    if (rank == 0)
        ALLDATA = MNIST17();

    timestamp_type start, stop;
    int t;
    int s;
    double e;
    double Z;
    double *error = malloc(T*sizeof(double));
    double *alpha = malloc(T*sizeof(double));
    double *running_error = malloc(T*sizeof(double));
    double sum;
    double send_weights[N];    

    Node **H = malloc(T*sizeof(Node*));
    for (t = 0; t < T; ++t) {
        H[t] = malloc(sizeof(Node));
        H[t]->parent = NULL;
    }

/*
    if (rank == 0) {
        printf("testing labels\n");
        for ( i = 0; i < n; ++i) {
            if (ALLDATA[i][D-1]*base[i]->label < 1)
                printf("WRONGLABEL i = %d", i);
        }
    }
*/
    
    printf("Starting AdaBoost\n");
    MPI_Barrier(MPI_COMM_WORLD);
    get_timestamp(&start);

    //int spaces;////////////////////////////////////////////////////

    for (t = 0; t < T; ++t) {

        //Pre-sort data before tree-building
        for (feat = 0; feat < num_features; ++feat) {
            PodSort(data[feat], 0, n-1, feat);
        }

        //Building tree
        if (rank == 0)
            printf("t = %d: Building tree\n", t); 
        ParallelSplit(H[t], data, n, 0, 0, rank, num_features);

        if (rank == 0) {
            e = PError(H[t], ALLDATA, base, n);
            error[t] = e;
            alpha[t] = 0.5*log((1 - e)/e);
            Z = 2*sqrt(e*(1 - e));
            for (i = 0; i < n; ++i) {
                base[i]->weight = base[i]->weight*exp(-alpha[t]*WeakLearner(H[t], ALLDATA[i])*base[i]->label)/Z;
                send_weights[i] = base[i]->weight;
            }

/*
            /////////////////////////////////////////////////////
            spaces = 0;
            for (i = 0; i < n; ++i) {
                if (WeakLearner(H[t], ALLDATA[i])*base[i]->label < 0) {
                    printf("%6d ", i);
                    spaces++;
                    if (spaces==14) {
                        printf("\n");
                        spaces = 0;
                    }
                }
            }
            ////////////////////////////////////////////////////
*/

            running_error[t] = 0;
            for (i = 0; i < n; ++i) {
                sum = 0;
                for (s = 0; s <= t; ++s)
                    sum += alpha[s]*WeakLearner(H[s], ALLDATA[i]);

                if (base[i]->label*sum < 0)
                    running_error[t] += 1./n;
            }

            printf("error = %f\n", running_error[t]);
        }
        
        MPI_Bcast(send_weights, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (i = 0; i < n; ++i)
            base[i]->weight = send_weights[i];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    get_timestamp(&stop);
    
    double elapsed = timestamp_diff_in_seconds(start, stop);
    printf("Elapsed time: %f seconds\n", elapsed);

    free(error);
    free(alpha);
    for (t = 0; t < T; ++t)
        TreeFree(H[t]);
    free(H);
    free(running_error);

    for (i = 0; i < n; ++i)
        free(base[i]);
    free(base);

    for (feat = 0; feat < num_features; ++feat)
        free(data[feat]);
    free(data);

    if (rank == 0) {
        for (i = 0; i < n; ++i)
            free(ALLDATA[i]);
    }
    free(ALLDATA);

    MPI_Finalize();

    return 0;
}




//*/











