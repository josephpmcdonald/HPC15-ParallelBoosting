#include <stdio.h>
#include <stdlib.h>
#include <math.h>


const int classes = 2;
const int N = 1000;
const int d = 10;

struct node {
    struct node *p;
    int ind;
    double split;
    struct node *right;
    struct node *left;
};


int sort(float data[][d], int n, int a) {
    //Sort data on index a and return
    //
    return 0;

}


int split(float data[][d], int n, int a, int pos) {

/* Returns the index of the table with best threshold split.
 * Partition rows up to and including that index from
 * everything afterwards.
 *
 * data = node data
 * n = length of data
 * a = index to split on
 * pos = number of positive labels
 *
 * Thus:
 * (pos - lpos) = right pos count
 * i = total left count
 * n - i = total right count
 *
 */

    int neg = n - pos;
    int lpos = 0;
    //int lneg = 0;
    double lp;
    double rp;

    double threshold;
    double P;
    double argmin;
    double threshmin;
    double Pmin;
    Pmin = 2.*(pos/n)*(neg/n); //initial purity of node

    //Tabulate purity for each possible threshold split    
    i = 0;
    while (i < n) {

        threshold = data[i][a];

        while (data[i][a] == threshold) {
            if (data[i][d-1] == 1)
                lpos += 1;
            ++i;
        }

        //right = n - i;
        //lneg = i - lpos;

        //Note that points on left = i, right = n-i
        lp = (double) lpos/i;
        rp = (double) (pos - lpos)/(n - i);

        //P = 2mu_l*lp(1-lp) + 2mu_r*rp(1-rp)
        P = 2.*lpos/n*(1 - lp) + 2.*(pos - lpos)/n*(1 - rp);

        //Split data at threshold/index with min purity
        if (P < Pmin) {
            Pmin = P;
            argmin = i - 1;
            threshmin = threshold;
        }
    }

    return argmin;
}


struct node *buildTree( ) {


    //Get the initial counts for positive/negatives.
    int i = 0;
    for (i = 0; i < n; ++i) {
        if (data[d*i + d-1] == 1.)
            pos += 1;
    }

    if (Pmin < 0.05 || n < 10) {
        struct node *leaf;
        

        return leaf;
    }

    int a;
    double Pmin = 100;
    double threshmin; 

    for (i = 0; i < d-1; ++a) {
        P = split(data, n, i, pos);
        if (P < Pmin) {
            a = i;
            Pmin = P;
            threshmin = data[a][];
        }
    }


    struct node *Node = malloc(sizeof *Node);
    Node->p = parent;
    Node->ind = a;
    Node->split = threshmin;

    //copy the data into     
    double *ldata = malloc(argmin * sizeof(double*));
    double *rdata = malloc((n - argmin) * sizeof(double*));
    double *lmatrix = malloc(argmin * d * sizeof(double));
    double rmatrix = malloc((n - argmin) * d * sizeof(double));
    for (i = 0; i < argmin; ++i)
        ldata[i] = &lmatrix[i*d];
    for (i = 0; i < n-argmin; ++i)
        rdata[i] = &rmatrix[i*d];

    int j = 0;
    for (i = 0; i < argmin; ++i) {
        for (j = 0; j < d; ++j)
            ldata[i][j] = data[i][j];
    }

    for (i = 0; i < n - argmin; ++i) {
        for (j = 0; j < d; ++j)
            rdata[i][j] = data[i + argmin][j];
    }
}






struct node* countsplit(float data[][d], int n, int a, node *) {

/* data = node data
 * n = length of data
 * a = index to split on
 *
 * Table computes/stores impurity for splitting with data <= threshold
 * counts[i][0] = ith threshold
 * counts[i][1] = left pos countcounts[i][2] = left neg count = i+1 - left pos count
 * counts[i][3] = impurity (Gini: 2*p*(1-p), entropy: p*log p + (1-p)log(1-p))
 *
 * Thus:
 * (pos - counts[i][1]) = right pos count
 * (neg - counts[i][2]) = right neg count
 * i+1 = total left count
 * n-(i+1) = total right count
 *
 */
    
    //
    //Counts stores pointer to each row of matrix
    double *matrix = calloc(n*3, sizeof(double));
    if (matrix)
        double *counts[] = malloc(n*sizeof(double*));

    int pos = 0;
    int neg = 0;
    int lpos = 0;
    int left = 0;
    int right = n;
    double lp;
    double rp;
    double purity;
    double argmin;
    double threshmin;

    //Initialize pointers in counts
    int i = 0;
    for (i = 0; i < n; ++i)
        counts[i] = &matrix[3*i];

    //Get the initial counts for positive/negatives.
    for (i = 0; i < n; ++i) {
        if (data[d*i + d-1] == 1.)
            pos += 1;
    }

    neg = n - pos;
    
    i = 0;
    int ind = 0;
    double threshold;
    double Pmin = 100.; //Upper bound on purity
    int argmin;
    while (i < n) {

        threshold = data[i][a];
        counts[ind][0] = threshold;

        while (data[i][a] == threshold) {
            if (data[i][d-1] == 1) {
                counts[ind][1] += 1;
                lpos += 1;
            }

            ++i;
            //Now left = i

            //else lneg = i - lpos

//            else if (data[i][d-1] == -1)
//                counts[ind][2] += 1;
//            else
//                printf("Data not clean.\n");

        }

        //Note that points on left = i, right = n-i
        lp = counts[ind][1]/i;
        //lp = (double) lpos/i;
        rp = (double) (pos - counts[ind][1])/(n - i);
        //rp = (double) (pos - lpos)/(n - i);

        //P = (left/n)*2*lp*(1 - lp) + (right/n)*2*rp*(1 - rp) and counts[ind][1]/n = (left/n)*lp
        counts[ind][2] = 2*(counts[ind][1]/n)*(1 - lp) + 2*(n - i)/n*rp*(1 - rp);
        purity = 2.*lpos/n*(1 - lp) + 2.*(pos - lpos)/n*(1 - rp);

        if (purity < Pmin) {
            Pmin = purity;
            argmin = i; //Note this is the index of the first right point
            threshmin = threshold;
        }

        //++ind;
    }
    
    for (i = 0; i < argmin; ++i) {
        

    }


    return 0;
}






int main(int argc, char *argv[]) {
    

    double a[d];
    a[0] = 1.;
    printf("%f\n", a[0]);
    return 0;

}






