#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <tree.h>


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
    double y_hat;

    for (i = 0; i < n; ++i) {
        if (WeakLearner(tree, data[i])*data[i][D+1] < 0)
            wrong++;
    }

    return (double) wrong/n;
}

double AdaBoost(double **data, int n) {
    
    int i;
    int t;
    int T = 100;
    double *error = malloc(T*sizeof(double));
    double *alpha = malloc(T*sizeof(double));
    for (i = 0; i < n; ++i) 
        data[i][D+1] = 1./n;

    //allocate array of pointers to tree nodes

    Node *H[] = malloc(T*sizeof(Node));
    for (t = 0; t < T; ++t) {
        H[i] = 
        error[t] = Error(H[i], data, n);
    }

    
}



