#ifndef IMPRIMIR_H_INCLUDED
#define IMPRIMIR_H_INCLUDED

//Fun��o para imprimir matrizes arbitr�rias
void ImprimeIJ(double **M, int maxI, int maxJ);
//Fun��o para imprimir matrizes char arbitr�rias
void ImprimeChar(char **M, int maxI, int maxJ);
//Fun��o para imprimir matrizes INT arbitr�rias (flags)
void ImprimeFlag(int **M, int maxI, int maxJ);

void Imprime3DTudo(double ***M, int maxI, int maxJ, int maxK);
void Imprime3DUm(double ***M, int maxI, int maxJ, int k);

#endif
