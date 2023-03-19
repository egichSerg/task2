#define max(X, Y) ((X) < (Y) ? (Y) : (X))

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    int netSize = 128;
    double minError = 0.000001;
    int maxIterations = 1000000;
    char* end;
    if (argc != 4){
        printf("You must enter excatly 3 arguments:\n1. Grid size (one number)\n2. Minimal error\n3. Iterations amount\n");
        return -1;
    }
    else{
        netSize = strtol(argv[1], &end, 10);
        minError = strtod(argv[2], &end);
        maxIterations = strtol(argv[3], &end, 10);
    }

    printf("%d, %0.20g, %d\n", netSize, minError, maxIterations);


    //values of net edges
    double tl = 10, //top left
        tr = 20, //top right
        bl = 30, //bottom left
        br = 20; //bottom right

    double horizontalStepTop = (tl - tr) / (netSize - 1), 
        horizontalStepBottom = (bl - br) / (netSize - 1), 
        verticalStepLeft = (tl - bl) / (netSize - 1), 
        verticalStepRight = (tr - br) / (netSize - 1);

    double** thermalConductivityMatrix = (double**)malloc(sizeof(double*) * netSize);
    double** thermalConductivityMatrixMod = (double**)malloc(sizeof(double*) * netSize);

#pragma acc data create(thermalConductivityMatrix[0:netSize][0:netSize], thermalConductivityMatrixMod[0:netSize][0:netSize])
{
    //init matrix
    for (int i = 0; i < netSize; i++)
    {
        thermalConductivityMatrix[i] = (double*)calloc(netSize, sizeof(double));
        thermalConductivityMatrixMod[i] = (double*)calloc(netSize, sizeof(double));
        if (i == 0 || i == netSize -1) {
            double step = verticalStepRight;
            double initNumber = tr;
            if (i == 0) { step = verticalStepLeft; initNumber = tl; }
            for (int j = 0; j < netSize; j++) {
                thermalConductivityMatrix[i][j] = initNumber + step * j;
            }
            
        }
        else {
            thermalConductivityMatrix[i][0] = tl + horizontalStepTop * i; 
            thermalConductivityMatrix[i][netSize - 1] = bl + horizontalStepBottom * i;
            }
        thermalConductivityMatrixMod[i][0] = tl + horizontalStepTop * i;
        thermalConductivityMatrixMod[i][netSize - 1] = bl + horizontalStepBottom * i;
    }


    double error = 10.;
    int iteration = 0;

    while (error > minError && iteration < maxIterations) {
        error = 0.;

        #pragma acc parallel loop reduction(max:error) collapse(2) independent
        for (int i = 1; i < netSize - 1; i++) {
            for (int j = 1; j < netSize - 1; j++) {
                thermalConductivityMatrixMod[i][j] = 0.25 * (
                    thermalConductivityMatrix[i][j + 1] + 
                    thermalConductivityMatrix[i][j - 1] + 
                    thermalConductivityMatrix[i + 1][j] + 
                    thermalConductivityMatrix[i - 1][j]);
                error = max(error, thermalConductivityMatrixMod[i][j] - thermalConductivityMatrix[i][j]);
            }
            
        }
        
        #pragma acc parallel loop collapse(2) independent
        for (int i = 0; i < netSize; i++) {
            for (int j = 0; j < netSize; j++) {
                thermalConductivityMatrix[i][j] = thermalConductivityMatrixMod[i][j];
            }
            
        }      
        iteration++;

        if (iteration % 1000 == 0) printf("iteration: %d error = %0.20g\n", iteration, error);
    }
    printf("Final error: %0.20g", error);
}
    return 0;
}