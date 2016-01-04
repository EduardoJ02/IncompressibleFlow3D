#include <stdlib.h>
#include <math.h>
#include "Alocar.h"

//Função para alocar matrizes iM x jM
//Function to alloc iMxjM arrays
double **Matriz(int iM, int jM){
    double **M; int i;

    M = malloc(iM * sizeof(double *));
    for (i = 0; i < iM; i++)
        M[i] = malloc(jM * sizeof(double));

    return M;
}

//Função para alocar matrizes iM x jM x kM
//Function to alloc iMxjMxkM arrays
double ***Matriz3D(int iM, int jM, int kM){
    double ***M; int i, j;

    M = (double ***)malloc(iM * sizeof(double **));

    for (i = 0; i < iM; i++){
        M[i] = (double **) malloc(jM * sizeof(double*));

        for (j = 0; j < jM; j++){
            M[i][j] = (double *) malloc(kM * sizeof(double));
        }
    }

    return M;
}

//Função para alocar matriz inteira iM x jM para flags
//Function to alloc iM x jM int flag array
int **FlagMatriz(int iM, int jM){
    int **M; int i;

    M = malloc(iM * sizeof(int *));
    for (i = 0; i < iM; i++)
        M[i] = malloc(jM * sizeof(int));

    return M;
}


//Função para alocar matriz inteira iM x jM x kM para flags
//Function to alloc iM x jM x kM int flag array
int ***FlagMatriz3D(int iM, int jM, int kM){
    int ***M; int i, j;

    M = (int ***)malloc(iM * sizeof(int **));

    for (i = 0; i < iM; i++){
        M[i] = (int **) malloc(jM * sizeof(int*));

        for (j = 0; j < jM; j++){
            M[i][j] = (int *) malloc(kM * sizeof(int));
        }
    }

    return M;
}


//Função para alocar matriz iM x jM de char
//Function to alloc iM x jM char array
char **CharMatriz(int iM, int jM){
    char **M; int i;

    M = (char**) malloc(iM * sizeof(char *)*iM);
    for (i = 0; i < iM; i++)
        M[i] = (char*)malloc(jM * sizeof(char));

    return M;
}

//Função para calcular o elemento máximo de uma matriz
    double Maximo(double ***M, int iM, int jM, int k){
        int i, j;
        double Max=M[0][0][0];
        for (i = 0; i < iM; i++) {
              for (j = 0; j < jM; j++){
                  if (Max < fabs(M[i][j][k])){
                     Max = fabs(M[i][j][k]);
                  }
              }
        }
        return Max;
    }
