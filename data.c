#include <stdlib.h>
#include <stdio.h>
#include "tree.h"

/* Functions for data IO, both retrieval and allocation. We choose to store
 * data in arrays of doubles. If D represents the number of features for a
 * particular observation x, then x[0] through x[D-1] will store these
 * features, x[D] will store the label, and x[D+1] will store the weight
 * for the observation as used in the AdaBoost algorithm.
 *
 * For MNIST17 and MNIST49, images are 28x28 arrays.
 *
 * Counts: 1/7=13007, 4/9=11791, other=35202    
 */


double **MNIST17() {

    int i;
    int j;
    int n = 13007;
    unsigned char c;
    FILE *labels;
    labels = fopen("train-labels.idx1-ubyte", "rb");
    FILE *images;
    images = fopen("train-images.idx3-ubyte", "rb");

/*
    //Print first 30 characters and first 5 images, check identifiable screen pattern
    fseek(labels, 8, SEEK_SET);
    for (i = 0; i < 30; ++i) {
        fread(&c, 1, 1, labels);
        printf("%d\n", c);
    }
    fseek(images, 16+28*28*0, SEEK_SET);
    for (i = 0; i < 28*5; ++i) {
        for (j = 0; j < 28; ++j) {
            fread(&c, 1, 1, images);
            printf("%4d ", c);
        }
        printf("\n");
    }
*/

    //Point to beginning of labels and images (after header info)
    fseek(labels, 8, SEEK_SET);
    fseek(images, 16, SEEK_SET);
    int count17 = 0;
    unsigned char label;
    unsigned char pixel;

    //Allocate pointers and array for each observation
    double **data17 = malloc(n*sizeof(double*));

    for (i = 0; i < n; ++i)
        data17[i] = malloc((28*28+2)*sizeof(double));

    //Count the occurrences of 1/7, and record image and label
    for (i = 0; i < 60000; ++i) {

        fread(&label, 1, 1, labels);

        if (label == 1 || label == 7) {
            //jump to correct line of pixel file
            fseek(images, 16+28*28*i, SEEK_SET);
            for (j = 0; j < 28*28; ++j) {
                fread(&pixel, 1, 1, images);
                data17[count17][j] = (double) pixel;
            }
            //record label at end of entry
            data17[count17][28*28] = (double) label;

            ++count17;
        }
    }

    printf("1/7:%d\n", count17);

    long lsize;
    fseek(images, 0, SEEK_END);
    lsize=ftell(images);
    //printf("The file is %d bytes\n", lsize);

    return data17;
}


double **MNIST49() {

    int i;
    int j;
    int n = 11791;
    unsigned char c;
    FILE *labels;
    labels = fopen("train-labels.idx1-ubyte", "rb");
    FILE *images;
    images = fopen("train-images.idx3-ubyte", "rb");

/*
    //Print first 30 characters and first 5 images, check identifiable screen pattern
    fseek(labels, 8, SEEK_SET);
    for (i = 0; i < 30; ++i) {
        fread(&c, 1, 1, labels);
        printf("%d\n", c);
    }
    fseek(images, 16+28*28*0, SEEK_SET);
    for (i = 0; i < 28*5; ++i) {
        for (j = 0; j < 28; ++j) {
            fread(&c, 1, 1, images);
            printf("%4d ", c);
        }
        printf("\n");
    }
*/

    //Point to beginning of labels and images (after header info)
    fseek(labels, 8, SEEK_SET);
    fseek(images, 16, SEEK_SET);
    int count49 = 0;
    unsigned char label;
    unsigned char pixel;

    //Allocate pointers and array for each observation
    double **data49 = malloc(n*sizeof(double*));

    for (i = 0; i < n; ++i)
        data49[i] = malloc((28*28+2)*sizeof(double));

    //Count the occurrences of 1/7, and record image and label
    for (i = 0; i < 60000; ++i) {

        fread(&label, 1, 1, labels);

        if (label == 4 || label == 9) {
            //jump to correct line of pixel file
            fseek(images, 16+28*28*i, SEEK_SET);
            for (j = 0; j < 28*28; ++j) {
                fread(&pixel, 1, 1, images);
                data49[count49][j] = (double) pixel;
            }
            //record label at end of entry
            data49[count49][28*28] = (double) label;

            ++count49;
        }
    }

    printf("4/9:%d\n", count49);

    long lsize;
    fseek(images, 0, SEEK_END);
    lsize=ftell(images);
    //printf("The file is %d bytes\n", lsize);

    return data49;
}



