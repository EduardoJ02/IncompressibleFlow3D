#include <stdlib.h>
#include <math.h>
#include "computaP.h"


//Gauss-Seidel method to solve the pressure equation
//Método de Gauss-Seidel para resolver a equação de pressão
void computaP(double ***mP3D, double ***mdP3D, double ***mP_ant3D, double ***residual3D, double ***mF3D, double ***mG3D, double ***mH3D, double ***u3D, double ***v3D, double ***w3D, int k){

    int it, i, j, cont=0;
    double a, b, c, d, ia, rMax, PSupp, PInfp, PDirp, PEsqp, PFrentep, PAtrasp;

    rMax=1;

    a=2*dt*(idx2+idy2+idz2);
    ia=1/a;

//    for (it=1; it<2; it++){
do{
        for (i = 0; i < iM; i++) {
          for (j = 0; j < jM; j++)
             mP_ant3D[i][j][k] = mP3D[i][j][k];
        }

            for (i = 0; i < iM; i++) {
                for (j = 0; j < jM; j++){

                    if (flag3D[i][j][k]==0){

                            b=(-dt*idx2)*(PSup + PInf);
                            c=(-dt*idy2)*(PDir + PEsq);
                            d=(-dt*idz2)*(PFrente + PAtras);
                            mP3D[i][j][k]=-ia*(b+c+d+mdP3D[i][j][k]);
                    }else{

                        if(Esquerdas){
                            PEsqp=PEsqC;
                        }else{
                            PEsqp=PEsq;
                        }
                        //printf("\nPEsqp=%.4f", PEsqp);


                        if(Direitas){
                            PDirp=PDirC;
                        }else{
                            PDirp=PDir;
                        }
                        //printf(" PDirp=%.4f", PDirp);

                        if(Inferiores){
                            PInfp=PInfC;
                        }else{
                            PInfp=PInf;
                        }
                        //printf(" PInfp=%.4f", PInfp);

                        if(Superiores){
                            PSupp=PSupC;
                        }else{
                            PSupp=PSup;
                        }
                        //printf(" PSupp=%.4f", PSupp);

                        if(Frentes){
                            PFrentep=PFrenteC;
                        }else{
                            PFrentep=PFrente;
                        }
                        //printf(" PFrentep=%.4f", PFrentep);

                        if(Atrases){
                            PAtrasp=PAtrasC;
                        }else{
                            PAtrasp=PAtras;
                        }
                        //printf(" PAtrasp=%.4f", PAtrasp);

                        //printf("\nPDirp=%.4f PSupp=%.4f PInfp=%.4f PEsqp=%.4f PFrentep=%.4f PAtrasp=%.4f", PDirp, PSupp, PInfp, PEsqp, PFrentep, PAtrasp);
                        b=(-dt*idx2)*(PSupp + PInfp);
                        c=(-dt*idy2)*(PDirp + PEsqp);
                        d=(-dt*idz2)*(PFrentep + PAtrasp);
                        mP3D[i][j][k]=-ia*(b+c+d+mdP3D[i][j][k]);

                        //printf("mP3D[%d][%d][%d]: %.4f\n", i, j, k, mP3D[i][j][k]);
                    }
                }
            }
            //printf("Matriz P depois de iteracao:\n");
            //Imprime3DUm(mP3D, iM,jM,k);

        rMax=0;
        for (i = 0; i < iM; i++) {
          for (j = 0; j < jM; j++){
             residual3D[i][j][k]=fabs(mP_ant3D[i][j][k]-mP3D[i][j][k]);
             if(residual3D[i][j][k]>rMax){
                rMax=residual3D[i][j][k];
             }


          }
        }
        //printf("rMax = %.4f\n", rMax);
        //printf("it = %d\n", cont+1);
        //cont++;

    }while (rMax>0.001);
    //}while (rMax>0.00004);
    //
    //}while (cont<4);

    //free(mP_ant);
}
