#include <stdlib.h>
#include <stdio.h>
#include "tree.h"

double **MNIST17() {
    //check that the data makes sense and matches labels, prints identifiable pattern on screen
    int i;
    int j;
    unsigned char c;
    FILE *labels;
    labels = fopen("train-labels.idx1-ubyte", "rb");
    fseek(labels, 8, SEEK_SET);
    for (i = 0; i < 30; ++i) {
        fread(&c, 1, 1, labels);
        printf("%d\n", c);
    }

    FILE *images;
    images = fopen("train-images.idx3-ubyte", "rb");
    fseek(images, 16+28*28*0, SEEK_SET);
    for (i = 0; i < 28*5; ++i) {
        for (j = 0; j < 28; ++j) {
            fread(&c, 1, 1, images);
            printf("%4d ", c);
        }
        printf("\n");
    }

    //Point to beginning of labels and images (after header info)
    fseek(labels, 8, SEEK_SET);
    fseek(images, 16, SEEK_SET);
    int count17 = 0;
    int count49 = 0;
    int count00 = 0;
    unsigned char label;
    unsigned char pixel;


    double **data17 = malloc(13007*sizeof(double*));
    double **data49 = malloc(11791*sizeof(double*));
    double **others = malloc(35202*sizeof(double*));

    for (i = 0; i < 13007; ++i)
        data17[i] = malloc((28*28+1)*sizeof(double));
    for (i = 0; i < 11791; ++i)
        data49[i] = malloc((28*28+1)*sizeof(double));
    for (i = 0; i < 35202; ++i)
        others[i] = malloc((28*28+1)*sizeof(double));

    //Count the occurrences of 1/7, 4/9, or others
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
        else if (label == 4 || label == 9) {

            //fseek(images, 16+28*28*i, SEEK_SET);
            for (j = 0; j < 28*28; ++j) {
                fread(&pixel, 1, 1, images);
                data49[count49][j] = (double) pixel;
            }
            data49[count49][28*28] = (double) label;

            ++count49;
        }
        else {

            //fseek(images, 16+28*28*i, SEEK_SET);
            for (j = 0; j < 28*28; ++j) {
                fread(&pixel, 1, 1, images);
                others[count00][j] = (double) pixel;
            }
            others[count00][28*28] = (double) label;

            ++count00;
        }
    }

    //1/7: 13007, 4/9: 11791, other: 35202    
    printf("1/7:%d\n", count17);
    printf("4/9:%d\n", count49);
    printf("other:%d\n", count00);

    /*
    printf("data49[0]=%f\n", data49[0][28*28]);
    printf("data49[1]=%f\n", data49[1][28*28]);
    printf("others[0]=%f\n", others[0][28*28]);
    printf("others[1]=%f\n", others[1][28*28]);
    */

    long lsize;
    fseek(images, 0, SEEK_END);
    lsize=ftell(images);
    //printf("The file is %d bytes\n", lsize);

    for (i = 0; i < 11791; ++i)
        free(data49[i]);
    free(data49);
    for (i = 0; i < 35202; ++i)
        free(others[i]);
    free(others);

    return data17;
}

int main(int argc, char *argv[]) {

    double **data17 = MNIST17();

    int i;
    for (i = 0; i < 10; ++i)
        printf("data17[%d]=%f\n", i, data17[i][28*28]);

    for (i = 0; i < 13007; ++i)
        free(data17[i]);
    free(data17);

    return 0;
}



