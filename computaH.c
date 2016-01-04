#include <stdlib.h>
#include "computaH.h"
#include "Discret.h"

/**
Stencil for H in a staggered grid
Estêncil para H em uma malha deslocada

k=3/2
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |       |   w   |       |       j
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


k=1
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|---v---|-------|-------
       |       |       |       |
       |       u       u       |       j
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
       |       |   w   |       |       j+1
       |       |       |       |
-------|-------|-------|-------|-------
       |       |       |       |
       |   w   |   w*  |   w   |       j
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


k=0
       |       |       |       |
       |       |       |       |       j+1
       |       |       |       |
-------|-------|---v*--|-------|-------
       |       |       |       |
       |       u       u*      |       j
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


double H3D(int i, int j, int k, double ***u3D, double ***v3D, double ***w3D){

    double Hf;

    Hf = w3D[i][j][k] + dt*(  iRe*laplaceW(i, j, Aai, Bb, Cci, Aaj, Ccj, Aa, Cc)
                              -uw_x(i,j, B, Bf, Bb, Aai, C, Cf, Cci)
                              -vw_y(i, j, b, bf, Bb, Aaj, e, ef, Ccj)
                              -w2(i,j, Bb, Aa, Cc) );
    return Hf;
}


double H3Dc(int i, int j, int k, double ***u3D, double ***v3D, double ***w3D){

    double Hf, Aaip, Bbp, Ccip, Aajp, Ccjp, Aap, Ccp, Bp, Bfp, Cp, Cfp, bp, bfp, ep, efp;


    //Sem problemas
    Bbp=Bb;
    Bp=B;
    bp=b;


    if (Frentes){/**k+1**/

        Aap=-Cc;
        bfp=-b;
        Bfp=-B;
    }else{
        Aap=Aa;
        bfp=bf;
        Bfp=Bf;
    }

    //Cc w3D[i][j][k-1]
    if(Atrases){/**k-1**/
        Ccp=-w3D[i][j][k+1];
    }else{
        Ccp=Cc;
    }


    //ef v3D[i][j-1][k+1]
    if (Frentes){/**k+1**/
        if ((Frentes) && (Inferiores)){/**k+1**//**j-1**/
            efp=v3D[i][j+1][k-1];
        }else{/**k+1**/
            efp=-v3D[i][j-1][k-1];
        }
    }else if (Inferiores){/**j-1**/
        efp=-v3D[i][j+1][k+1];
    }else{
        efp=ef;
    }

    //Cf u3D[i-1][j][k+1]
    if (Frentes){/**k+1**/
        if ((Frentes) && (Esquerdas)){/**k+1**//**i-1**/
            Cfp=u3D[i+1][j][k-1];
        }else{/**k+1**/
            Cfp=-u3D[i-1][j][k-1];
        }
    }else if (Esquerdas){/**i-1**/
        Cfp=-u3D[i+1][j][k+1];
    }else{
        Cfp=Cf;
    }

    //Aai w3D[i+1][j][k]
    if (Direitas){/**i+1**/
        Aaip=-Bb;
    }else{
        Aaip=Aai;
    }

    //Aaj w3D[i][j+1][k]
    if (Superiores){/**j+1**/
        Aajp=-Bb;
    }else{
        Aajp=Aaj;
    }

    //e v3D[i][j-1][k]
    //Ccj w3D[i][j-1][k]
    if (Inferiores){/**j-1**/
        ep=-v3D[i][j+1][k];
        Ccjp=-Bb;
    }else{
        ep=e;
        Ccjp=Ccj;
    }

    //C u3D[i-1][j][k]
    //Cci w3D[i-1][j][k]
    if (Esquerdas){/**i-1**/
        Cp=-u3D[i+1][j][k];
        Ccip=-Bb;
    }else{
        Cp=C;
        Ccip=Cci;
    }


    Hf = w3D[i][j][k] + dt*(  iRe*laplaceW(i, j, Aaip, Bbp, Ccip, Aajp, Ccjp, Aap, Ccp)
                              -uw_x(i,j, Bp, Bfp, Bbp, Aaip, Cp, Cfp, Ccip)
                              -vw_y(i, j, bp, bfp, Bbp, Aajp, ep, efp, Ccjp)
                              -w2(i,j, Bbp, Aap, Ccp) );
    return Hf;
}

