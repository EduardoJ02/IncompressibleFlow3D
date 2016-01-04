#ifndef ALOCAR_H_INCLUDED
#define ALOCAR_H_INCLUDED

//Fun��o para alocar matrizes iM x jM
double **Matriz(int iM, int jM);
//Fun��o para alocar matriz inteira para flags
int **FlagMatriz(int iM, int jM);
//Fun��o para alocar matriz de char
char **CharMatriz(int iM, int jM);

int ***FlagMatriz3D(int iM, int jM, int kM);
double ***Matriz3D(int iM, int jM, int kM);

double Maximo(double ***M, int iM, int jM, int k);

#endif
