#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tree.h"

//Note that D represents the number of features + 1 (for the label). 

double WeakLearner(Node *tree, double *data) {

/* WeakLearner is a function that takes an observation/vector of features and
 * returns a prediction on the label of that observation. The parameter tree is
 * a pointer to the root of the decision tree that WeakLearner represents as a
 * base classifier.
 *
 */

    double sign;
    

    return sign;
}


double Error(Node *tree, double **data, int n) {

    int i;
    int wrong = 0;

    for (i = 0; i < n; ++i) {
        if (WeakLearner(tree, data[i])*data[i][D-1] < 0)
            wrong++;
    }

    return (double) wrong/n;
}


double AdaBoost(double **data, int n) {
    
    int i;
    int t;
    int T = 1;
    double e;
    double Z;
    double final_error;
    double *error = malloc(T*sizeof(double));
    double *alpha = malloc(T*sizeof(double));
    Node **H = malloc(T*sizeof(Node*));
    for (t = 0; t < T; ++t)
        H[t] = malloc(sizeof(Node));

    for (i = 0; i < n; ++i) 
        data[i][D] = 1./n;

    for (t = 0; t < T; ++t) {
        BuildTree(H[t], data, n);
        e = Error(H[t], data, n);
        error[t] = e;
        alpha[t] = 0.5*log((1 - e)/e);
        Z = 2*sqrt(e*(1 - e));
        for (i = 0; i < n; ++i){
            data[i][D] = data[i][D]*exp(-alpha[t]*WeakLearner(&H[t], data[i])*data[i][D-1])/Z;
        }
    }

    double sum = 0;
    double y_hat;

    free(error);
    free(alpha);
    for (t = 0; t < T; ++t)
        TreeFree(H[t]);
    free(H);

    return final_error;
}



