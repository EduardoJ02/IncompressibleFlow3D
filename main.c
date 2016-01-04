#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Alocar.h"
#include "Imprimir.h"
#include "Arquivos.h"
#include "Discret.h"
#include "computaF.h"
#include "computaG.h"
#include "computaH.h"
#include "computaP.h"


double xlength, ylength, zlength, dx, dy, dz, idx2, idy2, idx, idy,idz, idz2, dt, idt;
int i,j, k, iM, jM, kM, n=0;
double Re, iRe, gama, eps, U_s;
int **flag;int ***flag3D;

int main(){

    double t=0, t_f=0.6, alfa, omega;
    int it, itM;

    //Lê as dimensões do problema a partir do arquivo
    //Reads the problem dimension from file
    leArquivoDim();
    printf("\n\nLinhas iM = %d\n", iM);
    printf("Colunas jM = %d", jM);
    kM=iM;


/**Variáveis do Problema****************************************************************************************/
    xlength=1.0; ylength=1.0; zlength=1.0;
    dt=0.02; alfa=0.25;
    eps=0.001; omega=1.7; gama=0.9;
    itM=100;
    Re=1000;
    U_s=1;

    dx=xlength/iM;
    dy=ylength/jM;
    dz=zlength/kM;

    idx2=1/(dx*dx);
    idy2=1/(dy*dy);
    idz2=1/(dz*dz);
    idx=1/dx;
    idy=1/dy;
    idz=1/dz;
    idt=1/dt;
    iRe=1/Re;
/**Variáveis do Problema****************************************************************************************/

/**Alocação das Matrizes Iniciais*******************************************************************/
    //Allocates initial arrays
    flag=FlagMatriz(iM, jM);
    flag3D=FlagMatriz3D(iM, jM, kM);

    leArquivo();

    for (k=0;k<kM;k++){
        for (i = 0; i <iM; i++) {
          for (j = 0; j <jM; j++){
                if (k==0){
                    flag3D[i][j][k]=flag[i][j]+16;
                } else if (k==kM-1){
                    flag3D[i][j][k]=flag[i][j]+32;
                }else{
                    flag3D[i][j][k]=flag[i][j];
                }
          }
        }
    }

    printf("\nMatriz de Flags 3D\n\n");
    for (k=0;k<kM;k++){

        if(k==3){
            printf(".\n.\n.\n.");
        }

        if (k==0 || k==kM-1 || k==1 || k==kM-2){
            printf("\nk=%d\n\n",k);
            for (i = 0; i <iM; i++) {
              for (j = 0; j <jM; j++){
                printf("%d ", flag3D[i][j][k]);
              }
               printf("\n");
            }
        }
    }

/***************U*****************//***************U*****************/
    double **u;    u=Matriz(iM+1, jM);
    double **mF;   mF=Matriz(iM+1, jM);

    double ***u3D; u3D=Matriz3D(iM+1, jM, kM);
    double ***mF3D;mF3D=Matriz3D(iM+1, jM, kM);

    for (i = 0; i < iM+1; i++) {
      for (j = 0; j < jM; j++){
            u[i][j] = 0;
            mF[i][j] = 0;
        }
    }

    for (k = 0; k < kM; k++) {
        for (i = 0; i < iM+1; i++) {
          for (j = 0; j < jM; j++){
                u3D[i][j][k] = 0;
                mF3D[i][j][k] = 0;
            }
        }
    }

/***************U*****************//***************U*****************/
/***************V****************//***************V****************/
    double **v;    v=Matriz(iM, jM+1);
    double **mG;   mG=Matriz(iM, jM+1);

    double ***v3D; v3D=Matriz3D(iM, jM+1, kM);
    double ***mG3D;mG3D=Matriz3D(iM, jM+1, kM);

    //Valores iniciais TESTE
    for (i = 0; i < iM; i++) {
      for (j = 0; j < jM+1; j++)
         v[i][j] = 0;
         mG[i][j] = 0;
    }

    for (k = 0; k < kM; k++) {
        for (i = 0; i < iM; i++) {
          for (j = 0; j < jM+1; j++){
                v3D[i][j][k] = 0;
                mG3D[i][j][k] = 0;
            }
        }
    }
/***************V****************//***************V****************/
/***************W****************//***************W****************/

    double ***w3D;    w3D=Matriz3D(iM, jM, kM+1);
    double ***mH3D;   mH3D=Matriz3D(iM, jM, kM+1);

    for (k = 0; k < kM+1; k++) {
        for (i = 0; i < iM; i++) {
          for (j = 0; j < jM; j++){
                w3D[i][j][k] = 0;
                mH3D[i][j][k] = 0;
            }
        }
    }
/***************W****************//***************W****************/
/***************P****************//***************P****************/

    double **mP;    mP=Matriz(iM, jM);
    double **mdP;   mdP=Matriz(iM, jM);

    double ***mP3D;    mP3D=Matriz3D(iM, jM, kM);
    double ***mdP3D;   mdP3D=Matriz3D(iM, jM, kM);

    //Valores iniciais TESTE
    for (i = 0; i < iM; i++) {
      for (j = 0; j < jM; j++)
         mP[i][j] = 0;
         mdP[i][j]=0;
    }

    for (k = 0; k < kM; k++) {
        for (i = 0; i < iM; i++) {
          for (j = 0; j < jM; j++){
                mP3D[i][j][k] = 0;
                mdP3D[i][j][k] = 0;
            }
        }
    }
/**Alocação das Matrizes Iniciais*******************************************************************/

/**Determina dt************************************************************************************/

//Determines the time step according to the stability condition
//Determina um passo dt que satisfaça a condição de estabilidade
void passoTempo(void){

    for(k=0;k<kM;k++){

        double uMax;        uMax=Maximo(u3D, iM+1, jM, k);

        double vMax;        vMax=Maximo(v3D, iM, jM+1, k);


        double wMax;        wMax=Maximo(w3D, iM, jM, k);

        double dtmin;
        if ((dx)/abs(uMax) < (dy)/fabs(vMax) && uMax!=0 && vMax!=0){
            dtmin = (dx)/fabs(uMax);
        }else{
            dtmin = (dy)/fabs(vMax);
        }

        if ((dz)/abs(wMax) < dtmin && wMax!=0)
            dtmin = (dz)/fabs(wMax);

        if ((Re/2)*1/(idx2 + idy2 +idz2) < dtmin)
            dtmin=(Re/2)*1/(idx2 + idy2 + idz2);

        dt=alfa*dtmin;
        }
}
/**Determina dt************************************************************************************/

/**Computa F, G e dP*******************************************************************************/
void computaFGdP(void){

    //Computes F, G and H velocities in the interior and in the boundary points
    //Computa velocidades F, G e H nos pontos de interior e de contorno
    for (i = 0; i<iM; i++) {
      for (j = 0; j<jM; j++){
            if (flag3D[i][j][k]==0){/**F e G e H no interior**/
                mF3D[i][j][k]=F3D(i, j, k, u3D, v3D);
                mG3D[i][j][k]=G3D(i, j, k, u3D, v3D);
                mH3D[i][j][k]=H3D(i, j, k, u3D, v3D, w3D);
            } else {                /**F e G e H no contorno**/
                mF3D[i][j][k]=F3Dc(i, j, k, u3D, v3D);
                mG3D[i][j][k]=G3Dc(i, j, k, u3D, v3D);
                mH3D[i][j][k]=H3Dc(i, j, k, u3D, v3D, w3D);
            }

            mF3D[i][j][k]=mF3D[i][j][k]+F3DcExtra(i,j,k,u3D,v3D,w3D);
            mG3D[i][j][k]=mG3D[i][j][k]+G3DcExtra(i,j,k,u3D,v3D,w3D);
      }
    }

    //Computes the right-hand side of the pressure equation
    //Calcula o lado direito da equação de pressão
    for (i = 0; i<iM; i++) {
      for (j = 0; j<jM; j++){
            mdP3D[i][j][k]=(mF3D[i+1][j][k] - mF3D[i][j][k])*idx
                         + (mG3D[i][j+1][k] - mG3D[i][j][k])*idy
                         + (mH3D[i][j][k+1] - mG3D[i][j][k])*idz;
      }
    }

 }
/**Computa F, G e dP*******************************************************************************/

/**Calcula u e v no nível de tempo seguinte*************************/
    //Computes U and V in the next time level
    void Passo(void){

    for (i = 0; i < iM-1; i++) {
      for (j = 0; j < jM; j++){
            u3D[i+1][j][k] = mF3D[i+1][j][k] - dt*idx*(mP3D[i+1][j][k]-mP3D[i][j][k]);
      }
    }

    for (i = 0; i < iM; i++) {
      for (j = 0; j < jM-1; j++){
            v3D[i][j+1][k] = mG3D[i][j+1][k] - dt*idy*(mP3D[i][j+1][k]-mP3D[i][j][k]);
      }
    }

    for (i = 0; i < iM; i++) {
      for (j = 0; j < jM; j++){
            if(k<kM-1){
                w3D[i][j][k+1] = mH3D[i][j][k+1] - dt*idy*(mP3D[i][j][k+1]-mP3D[i][j][k]);
            }
      }
    }
}
/**Calcula u e v no nível de tempo seguinte*************************/

/***************Algoritmo  ********************/

double ***mP_ant3D;
mP_ant3D=Matriz3D(iM, jM, kM);
double ***residual3D;
residual3D=Matriz3D(iM, jM, kM);

while (t<t_f){
//while (n<120){


    //Can also run on a fixed time step dt as long as it satisfies the stability condition
    //Também pode-se tomar um passo de tempo dt fixo contanto que satisfaça a condição de estabilidade
    passoTempo();
    t=t+dt;
    printf("\nTempo: %2f",t);
    printf("\n");

    for(k=0;k<kM;k++) {
    //    k=1;

        computaFGdP();

        computaP(mP3D, mdP3D, mP_ant3D, residual3D, mF3D, mG3D, mH3D, u3D, v3D, w3D, k);

        Passo();
    }

    n++;
}

    escreveArquivo(t,u3D,v3D,w3D);

    return 0;
}

