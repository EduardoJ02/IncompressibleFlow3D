#ifndef IMPRIMIR_H_INCLUDED
#define IMPRIMIR_H_INCLUDED

//Função para imprimir matrizes arbitrárias
void ImprimeIJ(double **M, int maxI, int maxJ);
//Função para imprimir matrizes char arbitrárias
void ImprimeChar(char **M, int maxI, int maxJ);
//Função para imprimir matrizes INT arbitrárias (flags)
void ImprimeFlag(int **M, int maxI, int maxJ);

void Imprime3DTudo(double ***M, int maxI, int maxJ, int maxK);
void Imprime3DUm(double ***M, int maxI, int maxJ, int k);

#endif
