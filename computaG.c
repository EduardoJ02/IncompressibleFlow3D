#include <stdlib.h>
#include "computaG.h"
#include "Discret.h"

/**
Stencil for G in a staggered grid
Estêncil para G em uma malha deslocada


k=1
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |       j
       |       |       |       |
-------|-------|---v---|-------|-------
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
       |       |   w   |       |       j
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |   w*  |       |       j-1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |       |       |
       |       |       |       |
              i-1      i      i+1


k=0
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|---v---|-------|-------
       |       |       |       |
       |       u       u       |        j
       |       |       |       |
-------|---v---|---v*--|---v---|-------
       |       |       |       |
       |       u       u*      |       j-1
       |       |       |       |
-------|-------|---v---|-------|-------
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
       |       |   w   |       |       j
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |   w   |       |       j-1
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
       |       |       |       |       j
       |       |       |       |
-------|-------|---v---|-------|-------
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


//Computes G in the interior points
//Computa G nos pontos interiores
double G3D(int i, int j, int k, double ***u3D, double ***v3D){
    double Gf;

    Gf = v3D[i][j][k] + dt*( (iRe)*((v2_x2(i,j, a, b, c) + v2_y2(i,j, d, b, e)))
                   - uv_x(i,j, B, D, b, a, C, Q, c) - v2(i,j, b, d, e)  );
    return Gf;
}


//Computes G in the boundary points
//Computa G nos pontos nos contornos
double G3Dc(int i, int j, int k, double ***u3D, double ***v3D){
    double Gf, Bp, Cp, Dp, Qp, ap, bp, cp, dp, ep;

    Bp=B;
    bp=b;

    //printf("\ni=%d j=%d k=%d\n", i,j,k);

    if(Esquerdas){/**i-1**/
        //i-1
            //c=-b
            //Q=-u(i+1,j+1)
            //C=-u(i+1,j)

        cp=-b;
        Cp=-u3D[i+1][j][k];

    }
    else{
        cp=c;
        Cp=C;
    }


    if (Direitas){/**i+1**/
        //Direita
            //i+1
            //a=-b;
        ap=-b;
    }else{
        ap=a;
    }

    if (Superiores){/**j+1**/
        //Superior
            //j+1
            //d=-e;
            //Q=-C;
            //D=-B;
            //tenho d=v(i,j+1);
            dp=d;
            Dp=-B;
    }else{
        dp=d;
        Dp=D;
    }

    if (Inferiores){/**j-1**/
        //j-1
            //e=-d;
            ep=-d;
    }else{
        ep=e;
    }


    /**j+1**/ /**i-1**/
    if ((Esquerdas) && (Superiores)){
        Qp=u3D[i+1][j][k];
    }else if ( ((Esquerdas) && (Inferiores)) || (Esquerdas)){
        Qp=-u3D[i+1][j+1][k];
    }else if ( ((Superiores)&& (Esquerdas)) || (Superiores)){
        Qp=-C;
    }
    else{
        Qp=Q;
    }


    Gf = v3D[i][j][k] + dt*( (iRe)*((v2_x2(i,j, ap, bp, cp) + v2_y2(i,j, dp, bp, ep)))
                     - uv_x(i,j, Bp, Dp, bp, ap, Cp, Qp, cp) - v2(i,j, bp, dp, ep)  );
    return Gf;
}

//Computes G in the boundary points in the third dimension
//Computa G nos pontos nos contornos na terceira dimensão
double G3DcExtra(int i, int j, int k, double ***u3D, double ***v3D,double ***w3D){
    double Gf, bp, bfrentep, batrasp, Bbp, Ccp, Aajp, Qqjp;

    bp=b;
    Bbp=Bb;

    if(flag3D[i][j][k]==0){
        bfrentep=bfrente;
        batrasp=batras;
        Aajp=Aaj;
        Ccp=Cc;
        Qqjp=Qqj;

    }

    if(Frentes){
        bfrentep=-b;
        batrasp=batras;

        Ccp=Cc;
        if((Frentes) && (Superiores)){
            Aajp=-w3D[i][j-1][k];
            Qqjp=-w3D[i][j-1][k-1];
        }else{
            Aajp=Aaj;
            Qqjp=Qqj;
        }
    }


    if(Atrases){
        bfrentep=bfrente;
        batrasp=-b;

        Ccp=-w3D[i][j][k+1];

        if((Atrases) && (Superiores)){
            Aajp=-w3D[i][j-1][k];
            Qqjp=w3D[i][j-1][k+1];
        }else{
            Aajp=Aaj;
            Qqjp=-w3D[i][j+1][k+1];
        }
    }



    Gf =dt*(  (iRe)*(v2_z2(i,j, bfrentep, bp, batrasp)) - vw_z(i,j,Bbp, Aajp, bp, bfrentep, Ccp, Qqjp, batrasp));

    //printf("\nExtra G[%d][%d][%d] somou %.6f\n", i,j,k,Gf);

    return Gf;
}
