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


double Error(Node *tree, double **data, int n) {

/* Note that Error returns the weighted error on the data.
 *
 */

    int i;
    double error = 0;

    for (i = 0; i < n; ++i) {
        if (WeakLearner(tree, data[i])*data[i][D-1] < 0)
            error += data[i][D];
    }

    return error;
}


int main (int argc, char *argv[]) {

/*
 */

    MPI_Init();
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int num_features = (D-1)/p;
    if (rank < D % p)
        num_features += 1;

    //////COPYDATA/////////////////
    double **ALLDATA = MNIST17();

    //Allocate pointers to pods 
    Pod **data = malloc(n*sizeof(Pod*));
    //Each process stores array of pods containing feature
    for (i = 0; i < n; ++i)
        data[i] = malloc(sizeof(Pod));

    for (i = 0; i < n; ++i) {
        data[i]->key = i; //key
        data[i]->val = ALLDATA[i][rank]; //feature
        data[i]->label = ALLDATA[i][D-1]; //label
        data[i]->weight = 1./n; //weight
    }

    free(ALLDATA);
    ///////////////////////////////


    timestamp_type start, stop;
    int i;
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
    for (t = 0; t < T; ++t)
        H[t] = malloc(sizeof(Node));

    printf("Starting AdaBoost\n");
    get_timestamp(&start);
    for (i = 0; i < n; ++i) 
        data[i][D] = 1./n;

    for (t = 0; t < T; ++t) {
        printf("t = %d: Building tree\n", t);
        BuildTree(H[t], data, n);
        e = Error(H[t], data, n);
        error[t] = e;
        alpha[t] = 0.5*log((1 - e)/e);
        Z = 2*sqrt(e*(1 - e));
        for (i = 0; i < n; ++i){
            data[i][D] = data[i][D]*exp(-alpha[t]*WeakLearner(H[t], data[i])*data[i][D-1])/Z;
        }

        running_error[t] = 0;
        for (i = 0; i < n; ++i) {
            sum = 0;
            for (s = 0; s <= t; ++s)
                sum += alpha[s]*WeakLearner(H[s], data[i]);

            if (data[i][D-1]*sum < 0)
                 running_error[t] += 1./n;
        }

        printf("error = %f\n", running_error[t]);
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
        free(data17[i]);
    free(data17);

    MPI_Finalize();

    return 0;
}





