#include <stdlib.h>
#include "computaF.h"
#include "Discret.h"

/**
Stencil for F in a staggered grid
Estêncil para F em uma malha deslocada


k=1
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       u       |       j
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |       j-1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |
       |       |       |       |
              i-1      i      i+1

k=1/2
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |   w*  |   w   |       j
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |       j-1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |
       |       |       |       |
              i-1      i      i+1


k=0
       |       |       |       |
       |       |       u       |       j+1
       |       |       |       |
-------|-------|---v*--|---v---|-------
       |       |       |       |
       |       u       u*      u       j
       |       |       |       |
-------|-------|---v---|---v---|-------
       |       |       |       |
       |       |       u       |       j-1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |
       |       |       |       |
              i-1      i      i+1

k=-1/2
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |   w   |   w   |       j
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |       j-1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |
       |       |       |       |
              i-1      i      i+1


k=-1
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       u       |       j
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |       j-1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |
       |       |       |       |
              i-1      i      i+1

u*=u[i][j][k]
v*=v[i][j][k]
w*=w[i][j][k]

This code is set up with the boundary conditions for the Lid-Driven Cavity problem
Este código está configurado com as condições de contorno para o problema do Escoamento em uma Cavidade

Use this scheme to interpolate boundary velocities using the proper points in a staggered grid
Use este esquema para interpolar velocidades no contorno usando os pontos adequados na malha deslocada

**/

//Computes F in the interior points
//Computa F nos pontos interiores
double F3D(int i, int j, int k, double ***u3D, double ***v3D){
    double Ff;

    Ff = u3D[i][j][k] + dt*(  (iRe)*((u2_x2(i,j, A, B, C) + u2_y2(i,j, D, B, E)))
                                       - u2(i,j, B, A, C)  - uv_y(i,j, b, a, B, D, e, q, E)  );
    return Ff;
}

//Computes F in the boundary points
//Computa F nos pontos nos contornos
double F3Dc(int i, int j, int k, double ***u3D, double ***v3D){
    double Ff, Ap, Bp, Cp, Dp, Ep, ap, bp, ep, qp;

    Bp=B;
    bp=b;

    if(Esquerdas){/**i-1**/
        //Esquerda
        //i-1 dá problema
        //C=-A
        Cp=-A;
        //Ff=3;
    }else{
        Cp=C;
    }

    if (Direitas){/**i+1**/
            //Direita
            //i+1 dá problema
            //A=-C
            //a=-b;
            //q=-e;

            //tenho A=u(i+1,j);

        Ap=A;
        ap=-b;
        //qp=-e;
        //Ff=5;
    }else{
        Ap=A;
        ap=a;
        //qp=q;
    }


    if (Superiores){/**j+1**/
        //Superior
        //j+1
        //D=2*U_s-B
        if (flag3D[i][j][k]==Superior+Esquerda || flag3D[i][j][k]==Superior+Esquerda+Atras || flag3D[i][j][k]==Superior+Esquerda+Frente ){
                Dp=-B;
        } else{
                Dp=2*U_s-B;
        }
        //Ff=2;
    }else{
        Dp=D;
    }

    if (Inferiores){/**j-1**/
            //Inferior
            //j-1 dá problema
            //E=-B
            //e=-v(i,j+1);
            //q=-v(i+1,j+1);

        Ep=-B;
        ep=-v3D[i][j+1][k];
        //qp=-v3D[i+1][j][k];
        //Ff=1;
    }else{
        Ep=E;
        ep=e;
        //qp=q;
    }


    /**j-1**/ /**i+1**/
    if ((Direitas) && (Inferiores)){
        qp=v3D[i][j+1][k];
    }else if ( ((Esquerdas) && (Inferiores)) || (Inferiores) ){
        qp=-v3D[i+1][j+1][k];
    }else if ( ((Superiores)&& (Direitas)) || (Direitas)){
        qp=-e;
    }
    else{
        qp=q;
    }

    Ff = u3D[i][j][k] + dt*(  (iRe)*((u2_x2(i,j, Ap, Bp, Cp) + u2_y2(i,j, Dp, Bp, Ep)))
                                       - u2(i,j, Bp, Ap, Cp)  - uv_y(i,j, bp, ap, Bp, Dp, ep, qp, Ep)  );

    return Ff;
}

//Computes F in the boundary points in the third dimension
//Computa F nos pontos nos contornos na terceira dimensão
double F3DcExtra(int i, int j, int k, double ***u3D, double ***v3D, double ***w3D){
    double Ff, Bp, Bfrentep, Batrasp, Bbp, Ccp, Aaip, Qqip;

    Bp=B;
    Bbp=Bb;

    if(flag3D[i][j][k]==0){

        Bfrentep=Bfrente;
        Batrasp=Batras;
        Aaip=Aai;
        Ccp=Cc;
        Qqip=Qqi;
    }


    if(Frentes){
        //Bfrentep=-Batras;
        Bfrentep=-B;
        Batrasp=Batras;

        Ccp=Cc;
        if((Frentes) && (Direitas)){
            Aaip=-w3D[i-1][j][k];
            //Qqip=-w3D[i-1][j][k-1];
            Qqip=-w3D[i-1][j][k-1];
        }else{
            Aaip=Aai;
            Qqip=Qqi;
        }
    }

    if(Atrases){

        Bfrentep=Bfrente;
        Batrasp=-B;

        Ccp=-w3D[i][j][k+1];

        if((Atrases) && (Direitas)){
            Aaip=-w3D[i-1][j][k];
            Qqip=w3D[i-1][j][k+1];
        }else{
            Aaip=Aai;
            Qqip=-w3D[i+1][j][k+1];
        }

    }

    Ff =dt*(  (iRe)*(u2_z2(i,j, Bfrentep, Bp, Batrasp)) - uw_z(i,j,Bbp, Aaip, Bp, Bfrentep, Ccp, Qqip, Batrasp));

    //printf("\nExtra F[%d][%d][%d] somou %.6f\n", i,j,k,Ff);

    return Ff;
}
