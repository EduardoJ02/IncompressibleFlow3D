#include <stdlib.h>
#include <stdio.h>
#include "Arquivos.h"

#define Arquivo "Dom48x48.txt"

void Marcar(int **flag, char**flagchar,int iM, int jM){
    int i, j;
        for (i=0; i<=iM; i++) {
              for (j = 0; j<=jM; j++){

               if (flagchar[i][j]=='w'){
                    flag[j][iM-i]=3;
               }
               if (flagchar[i][j]=='s'){
                    flag[j][iM-i]=2;
               }
               if (flagchar[i][j]=='x'){
                    flag[j][iM-i]=10;
               }
               if (flagchar[i][j]=='d'){
                    flag[j][iM-i]=8;
               }
               if (flagchar[i][j]=='e'){
                    flag[j][iM-i]=1;
               }
               if (flagchar[i][j]=='y'){
                    flag[j][iM-i]=5;
               }
               if (flagchar[i][j]=='i'){
                    flag[j][iM-i]=4;
               }
               if (flagchar[i][j]=='z'){
                    flag[j][iM-i]=12;
               }
               if (flagchar[i][j]=='.'){
                    flag[j][iM-i]=0;
               }
               //fprintf(stderr, "%c", M[i][j]);
               //printf("\n");
            }
        }
    }

char** leArquivoDim(void){

        FILE *geometria;
        geometria = fopen(Arquivo, "r");

        int interioresCol=0, interioresLin=0, contt=0;

        int ch, n=0, cont=0,tam, i, j, k;

        if (geometria == NULL){
            printf("Erro\n");
            return NULL;
        }

        while ((ch = fgetc(geometria) != EOF)){
           cont++;
        }

        tam = cont;
        char *flagfile = (char*) calloc(tam, sizeof(char));
        rewind(geometria);

        printf("Geometria lida:\n\n");

        for(i=0; i<tam; i++)
           fscanf(geometria, "%c ", &flagfile[i]);

        j=0;i=0;
        for(k=0; k<tam; k++){

            if (contt==0 && flagfile[k]!=' ' && flagfile[k]!='l'){
                    interioresCol=interioresCol+1;
            }

           if (flagfile[k]=='l'){
                //printf("interiores = %d", interioresCol);
                interioresLin=interioresLin+1;

                contt=1;

                j=0;
                i=i+1;
                printf("\n");
           }
           else{

               fprintf(stderr, "%c", flagfile[k]);

               //flagchar[i+1][j+1]=flagfile[k];
               j=j+1;
           }
           iM=interioresLin;
           jM=interioresCol;


/*
            if (flagfile[k]=='f'){
                //printf("\n\nfim\n lidos %d caracteres do arquivo\n\n",k);
                printf("\nColunas iM = %d", iM);
                printf("\nLinhas jM = %d", jM);
                break;
           }
*/
        }

        fclose(geometria);

    }


char** leArquivo(void){

        FILE *geometria;
        geometria = fopen(Arquivo, "r");
        int i, j;

        //Cria matriz flag
        char **flagchar;
        flagchar=CharMatriz(iM+1,jM+1);

        //Valores iniciais TESTE
        for (i = 0; i <= iM; i++) {
          for (j = 0; j <= jM; j++)
                flagchar[i][j] = 'o';
        }


        int interioresCol=0, interioresLin=0;

        int ch, n=0, cont=0,tam, k;

        if (geometria == NULL){
            printf("Erro\n");
            return NULL;
        }

        while ((ch = fgetc(geometria) != EOF)){
           cont++;
        }

        tam = cont;
        char *flagfile = (char*) calloc(tam, sizeof(char));
        rewind(geometria);

        //printf("%d\n\n", tam);

        for(i=0; i<tam; i++)
           fscanf(geometria, "%c ", &flagfile[i]);

        j=0;
        i=0;
        for(k=0; k<tam; k++){

            if (flagfile[k]!=' ' && flagfile[k]!='l'){
                    //interioresCol=interioresCol+1;
                    flagchar[i][j]=flagfile[k];
                    //printf("ei: %c", flagchar[i][j]);


                    if (flagfile[k+1]=='l'){
                        flagchar[i][j]=flagfile[k];
                        j=0;
                        i=i+1;
                    }

                    j=j+1;
            }
            else{
                    j=0;
            }


            if (flagfile[k]=='f'){
                printf("\n\nfim\n lidos %d caracteres do arquivo\n\n",k);
                break;
           }
        }

        printf("\nMatriz de Flags inteira\n\n");
        Marcar(flag,flagchar,iM-1,jM-1);
        ImprimeFlag(flag,iM-1,jM-1);

        fclose(geometria);
    }


char** escreveArquivo(double t, double ***u3D, double ***v3D, double ***w3D){

        int i, j, k;

        FILE *saida;
        saida = fopen("SAIDA.txt", "w");

        FILE *saidaU;
        saidaU = fopen("U_3D.txt", "w");
        FILE *saidaV;
        saidaV = fopen("V_3D.txt", "w");
        FILE *saidaW;
        saidaW = fopen("W_3D.txt", "w");

        if (saida == NULL || saidaU == NULL || saidaV == NULL || saidaW == NULL){
            printf("Erro\n");
            return NULL;
        }

        fprintf(saidaU,"Geometria %d x %d x %d \nTempo = %.4fs \n\n", iM, jM, kM, t);
        fprintf(saidaV,"Geometria %d x %d x %d \nTempo = %.4fs \n\n", iM, jM, kM, t);
        fprintf(saidaW,"Geometria %d x %d x %d \nTempo = %.4fs \n\n", iM, jM, kM, t);



        for (k=0; k<kM; k=k+1){
            fprintf(saidaU,"U%d=[",k);
            for (i = 0; i <iM+1; i=i+1) {
              for (j = 0; j <jM; j=j+1){
                   if(j==0){
                    fprintf(saidaU, "0 ");
                   }

                    fprintf(saidaU, "%.8lf ", u3D[i][j][k]);
                  }
                  fprintf(saidaU," 1;\n");
              }

            fprintf(saidaU,"]\n\n");
        }



        for (k=0; k<kM; k=k+1){
            fprintf(saidaV,"V%d=[",k);
            for (i = 0; i<iM+1; i=i+1) {
              for (j = 0; j<jM+1; j=j+1){
                if(i==0){
                    fprintf(saidaV, "0 ");
                }
                else{
                    fprintf(saidaV, "%.8lf ", v3D[i-1][j][k]);
                }

                }
                fprintf(saidaV,";\n");
                }
             for (j = 0; j<jM+1; j=j+1){
                fprintf(saidaV, "0 ");
             }
             fprintf(saidaV,";\n");
             fprintf(saidaV,"]\n\n");
        }




        for (k=0; k<kM+1; k=k+1){
            fprintf(saidaW,"W%d=[",k);
            for (i = 0; i<iM; i=i+1) {
              for (j = 0; j<jM; j=j+1){
                /*if(i==0){
                    fprintf(saidaW, "0 ");
                }
                else{*/
                    fprintf(saidaW, "%.8lf ", w3D[i][j][k]);
                //}

                }
                fprintf(saidaW,";\n");
                }
             /*for (j = 0; j<jM+1; j=j+1){
                fprintf(saidaW, "0 ");
             }*/
             fprintf(saidaW,";\n");
             fprintf(saidaW,"]\n\n");
        }


/*
         fprintf(saida,"P=[");
        for (i = 0; i<iM; i=i+1) {
          for (j = 0; j<jM; j=j+1){
                fprintf(saida, "%.8lf ", mP[i][j]);
          }
          fprintf(saida,";\n");
        }
                fprintf(saida,";\n");
         fprintf(saida,"]\n\n");

*/

        fprintf(saida, "\niM=%d;jM=%d;\nz= linspace(0,1,%d);",iM,jM,iM+2);
        fprintf(saida, "\nplot(z,U(floor(iM/2),:));");
        fprintf(saida, "\nplot(z,V(:,floor(iM/2)));");

        fprintf(saida, "\nU(:,jM+1)=(U(:,jM+1)+U(:,jM+2))/2; U(:,jM+2)=[];");
        fprintf(saida, "V(iM+1,:)=(V(iM+1,:)+V(iM+2,:))/2; V(iM+2,:)=[]");


        fprintf(saida, "\n\nx = linspace(0,1,%d+1); y = linspace(0,1,%d+1);",iM, jM);
        //fprintf(saida, "\n\[X,Y]=meshgrid(x,y);\nscf(12); clf(12);\nchamp(x,y,U,V',rect=[0, 0,1,1], arfact=1.5);");
        fprintf(saida, "\n\[X,Y]=meshgrid(x,y);\nscf(12); clf(12);\nchamp1(x,y,U,V,rect=[-0.01, -0.01,1.01,1.01], arfact=0.5)");



      // fprintf(saida, "Texto: %lf %s\n", u, str);
       fclose(saida);

    }
