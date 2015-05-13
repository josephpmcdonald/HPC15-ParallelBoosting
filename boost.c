#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tree.h"
#include "util.h"

//Note that D represents the number of features + 1 (for the label). 


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

    int i;
    int wrong = 0;

    for (i = 0; i < n; ++i) {
        if (WeakLearner(tree, data[i])*data[i][D-1] < 0)
            wrong++;
    }

    return (double) wrong/n;
}


double AdaBoost(double **data, int n) {
    
    timestamp_type start, stop;
    int i;
    int t;
    int T = 1;
    double e;
    double Z;
    double *error = malloc(T*sizeof(double));
    double *alpha = malloc(T*sizeof(double));
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
        printf("checking error\n");
        e = Error(H[t], data, n);
        error[t] = e;
        alpha[t] = 0.5*log((1 - e)/e);
        Z = 2*sqrt(e*(1 - e));
        printf("reweighting\n");
        for (i = 0; i < n; ++i){
            data[i][D] = data[i][D]*exp(-alpha[t]*WeakLearner(H[t], data[i])*data[i][D-1])/Z;
        }
    }
    get_timestamp(&stop);
    double elapsed = timestamp_diff_in_seconds(start, stop);
    printf("Elapsed time: %f seconds\n", elapsed);

    double final_error = 0;
    double sum;
    for (i = 0; i < n; ++i) {
        sum = 0;
        for (t = 0; t < T; ++t)
            sum += alpha[t]*WeakLearner(H[t], data[i]);
        
        if (data[i][D-1]*sum < 0)
            final_error += 1./n;
    }


    free(error);
    free(alpha);
    for (t = 0; t < T; ++t)
        TreeFree(H[t]);
    free(H);

    return final_error;
}


int main(int argc, char *argv[]) { 


    double **data17 = MNIST17();
    //Node *root = malloc(sizeof(Node));
    //SplitNode(root, data17, 1000, 0, 0);

    //To check appearance of data17
    int i;
    int j;
    int k;
    for (i = 0; i < 10; ++i)
        printf("data17[%d]=%f\n", i, data17[i][D-1]);

    for (i = 0; i < 10; ++i) {
        for (j = 0; j < 28; ++j) {
            for (k = 0; k < 28; ++k)
                printf("%5.0f ", data17[i][28*j+k]);

            printf("\n");
        }

        printf("\n");
    }

    double error = AdaBoost(data17, 13007);
    printf("Final error: %f\n", error);

    for (i = 0; i < 13007; ++i)
        free(data17[i]);
    free(data17);



/*
    double DATA[20][D+1] = {{3, 6, 1, 0.05},
        {4, 1, -1, 0.05},
        {3, 1, -1, 0.05},
        {5, 7, 1, 0.05},
        {1, 6, 1, 0.05},
        {8, 2, -1, 0.05},
        {2, 8, 1, 0.05},
        {5, 2, -1, 0.05},
        {9, 1, -1, 0.05},
        {5, 4, -1, 0.05},
        {2, 1, -1, 0.05},
        {1, 6, 1, 0.05},
        {3, 9, 1, 0.05},
        {2, 1, -1, 0.05},
        {8, 3, -1, 0.05},
        {5, 1, -1, 0.05},
        {6, 4, -1, 0.05},
        {9, 5, 1, 0.05},
        {6, 0, -1, 0.05},
        {11, 11, 1, 0.05}};

    double *data[20];
    int i;
    int j;
    for (i = 0; i < 20; ++i)
        data[i] = &DATA[i][0];




    Sort(data, 0, 19, 1);

    double impurity;
    i = BestSplit(data, 20, 0, 1, 8, &impurity);
    printf("i = %d\n", i);


    Node *root = malloc(sizeof(Node));
    root->right = NULL;
    root->left = NULL;
    root->parent = NULL;
    BuildTree(root, data, 20);

    for (i = 0; i < 20; ++i) {
        for (j = 0; j < D; ++j)
            printf("%f ", data[i][j]);
        printf("\n");
    }

    TreePrint(root, 0);

    double pt[] = {2, 4, -1};
    double *point = &pt[0];
    printf("Test: f(pt) = %f\n", TestPoint(root, point));

    TreeFree(root);
*/

    return 0;
}




