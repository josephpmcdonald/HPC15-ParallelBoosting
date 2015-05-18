#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tree.h"
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
    int n = 13007;
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
    double **MYDATA = ParMNIST17(feature_list, num_features);

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

    ///////PRE-SORT///////////////////////////
    for (feat = 0; feat < num_features; ++feat)
        PodSort(data[feat], 0, n-1, feat);
    //////////////////////////////////////////

    free(MYDATA);

    double **ALLDATA = NULL;
    if (rank == 0)
        ALLDATA = MNIST17();

    /////////////////////////////////////////////////////////

    timestamp_type start, stop;
    int t;
    int s;
    int T = 50;
    double e;
    double Z;
    double *error = malloc(T*sizeof(double));
    double *alpha = malloc(T*sizeof(double));
    double *running_error = malloc(T*sizeof(double));
    double sum;

    Node **H = malloc(T*sizeof(Node*));
    for (t = 0; t < T; ++t) {
        H[t] = malloc(sizeof(Node));
        H[t]->parent = NULL;
    }


    printf("Starting AdaBoost\n");
    get_timestamp(&start);

    for (t = 0; t < T; ++t) {

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
            }

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
    }

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





