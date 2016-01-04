#include <stdlib.h>
#include "Imprimir.h"

    //Função para imprimir matrizes arbitrárias
    void ImprimeIJ(double **M, int maxI, int maxJ){
        int i=0;int j=0;

        for (i = 0; i <maxI; i++) {
          for (j = 0; j < maxJ; j++)
            printf("%.8lf ", M[i][j]);
           printf("\n");
        }
    }

    //Função para imprimir matrizes char arbitrárias
    void ImprimeChar(char **M, int maxI, int maxJ){
        int i; int j;
        for (i = 0; i <maxI; i++) {
          for (j = 0; j < maxJ; j++)
            printf("%c", M[i][j]);
           printf("\n");
        }
    }

    //Função para imprimir matrizes INT arbitrárias (flags)
    void ImprimeFlag(int **M, int maxI, int maxJ){
        int i; int j;
        for (i = 0; i <=maxI; i++) {
          for (j = 0; j <= maxJ; j++)
            printf("%d ", M[i][j]);
           printf("\n");
        }
    }

    void Imprime3DUm(double ***M, int maxI, int maxJ, int k){

        int i, j;
            printf("\nk=%d\n\n",k);
            for (i = 0; i <maxI; i++) {
              for (j = 0; j <maxJ; j++){
                printf("%.8f ", M[i][j][k]);
              }
               printf("\n");
            }
    }

    void Imprime3DTudo(double ***M, int maxI, int maxJ, int maxK){

        int i, j, k;

        for (k=0;k<maxK;k++){
            printf("k=%d\n\n",k);

            for (i = 0; i <maxI; i++) {
              for (j = 0; j <maxJ; j++){
                printf("%.8f ", M[i][j][k]);
              }
               printf("\n");
            }
        }
    }
