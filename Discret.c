#include <stdlib.h>
#include <math.h>
#include "Discret.h"

extern double dx, dy, idx2, idx, idy2, idy, idz, idz2, gama;

/**Discretizações para F**/
    double u2_x2(int i, int j, double Ap, double Bp, double Cp){
        double res;

        res=(Ap - 2.0*Bp + Cp)*idx2;

        return res;
    }

    double u2_y2(int i, int j, double Dp, double Bp, double Ep){
        double res;

        res=(Dp-2.0*Bp+Ep)*idy2;

        return res;
    }

    double u2(int i, int j, double Bp, double Ap, double Cp){
        double res;

        res=((((Bp + Ap))/2)*(((Bp + Ap))/2) - (((Bp + Cp))/2)*(((Bp+ Cp))/2))*idx    +
          gama*idx*0.25*( fabs(Bp+Ap)*(Bp-Ap)  - fabs(Cp+Bp)*(Cp-Bp));

        return res;
    }

    double uv_y(int i, int j, double bp, double ap, double Bp, double Dp, double ep, double qp, double Ep){
        double res;

        //res=((((bp + ap))/2)*(((Bp + Dp))/2) - (((ep + qp))/2)*(((Ep + Bp))/2))
         //+    gama*idy*0.25*(fabs(bp+ap)*(Bp-Dp) -  fabs(ep+qp)*(Ep-Bp));

        res=((bp + ap)*(Bp + Dp) +  gama*fabs(bp + ap)*(Bp - Dp)-
            (ep + qp)*(Ep + Bp) - gama*fabs(ep + qp)*(Ep - Bp))
       //             /(4.0*dy);
       *(0.25*idy); ///ARRUMEI

        return res;
    }

    //Para F Extra:
    double u2_z2(int i, int j, double Bfrentep, double Bp, double Batrasp){
        double res;

        res=(Bfrentep - 2.0*Bp + Batrasp)*idz2;

        return res;
    }

    double uw_z(int i, int j, double Bbp, double Aaip,double Bp,double Bfrentep,double Ccp,double Qqip,double Batrasp){
        double res;

        res=((Bbp + Aaip)*(Bp + Bfrentep) +  gama*fabs(Bbp + Aaip)*(Bp - Bfrentep) -
             (Ccp + Qqip)*(Batrasp + Bp)  -  gama*fabs(Ccp + Qqip)*(Batrasp - Bp))
             *(0.25*idz);

        return res;
    }


/**Discretizações para G**/
    double v2_x2(int i, int j, double ap, double bp, double cp){
        double res;

        res=(ap-2*bp+cp)*idx2;

        return res;
    }

    double v2_y2(int i, int j, double dp, double bp, double ep){
        double res;

        res=(dp-2*bp+ep)*idy2;

        return res;
    }

    double v2(int i, int j, double bp, double dp, double ep){
        double res;

        res=((((bp+ dp))/2)*(((bp+ dp))/2) - (((ep + bp))/2)*(((ep + bp))/2))*idy
        + gama*0.25*idy*(fabs(bp+dp)*(bp-dp)  -  fabs(ep+bp)*(ep-bp));

        return res;
    }

    double uv_x(int i, int j, double Bp, double Dp, double bp, double ap, double Cp, double Qp, double cp){
        double res;

         res = ((Bp + Dp)*(bp + ap)+
                    gama*fabs(Bp + Dp)*(bp - ap)-
                    (Cp + Qp)*(cp + bp)-
		    gama*fabs(Cp + Qp)*(cp - bp))
                    // /(4.0*dx);
                    *(0.25*idx); ///ARRUMEI

        return res;
    }

//Para G Extra:
    double v2_z2(int i, int j, double bfrentep, double bp, double batrasp){
        double res;

        res=(bfrentep - 2.0*bp + batrasp)*idz2;

        return res;
    }

    double vw_z(int i, int j, double Bbp, double Aajp,double bp,double bfrentep,double Ccp,double Qqjp,double batrasp){
        double res;

        res=((Bbp + Aajp)*(bp + bfrentep) +  gama*fabs(Bbp + Aajp)*(bp - bfrentep) -
             (Ccp + Qqjp)*(batrasp + bp)  -  gama*fabs(Ccp + Qqjp)*(batrasp - bp))
             *(0.25*idz);

        return res;
    }



/**Discretizações para H**/
double laplaceW(int i, int j, double Aaip,double Bbp,double Ccip,double Aajp,double Ccjp,double Aap,double Ccp){
        double res;

        res= (Aaip - 2.0*Bbp + Ccip)*idx2
            +(Aajp - 2.0*Bbp + Ccjp)*idy2
            +(Aap - 2.0*Bbp + Ccp)*idz2;

        return res;
}

double uw_x(int i, int j, double Bp,double Bfp,double Bbp,double Aaip,double Cp,double Cfp,double Ccip){
        double res;

        res=((Bp + Bfp)*(Bbp + Aaip) +  gama*fabs(Bp + Bfp)*(Bbp - Aaip) -
             (Cp + Cfp)*(Ccip + Bbp)  -  gama*fabs(Cp + Cfp)*(Ccip - Bbp))
             *(0.25*idx);

        return res;
}

double vw_y(int i, int j, double bp,double bfp,double Bbp,double Aajp,double ep,double efp,double Ccjp){
        double res;

        res=((bp + bfp)*(Bbp + Aajp) +  gama*fabs(bp + bfp)*(Bbp - Aajp) -
             (ep + efp)*(Ccjp + Bbp)  -  gama*fabs(ep + efp)*(Ccjp - Bbp))
             *(0.25*idy);

        return res;
}

double w2(int i, int j, double Bbp,double Aap,double Ccp){
        double res;

        res=((((Bbp+ Aap))/2)*(((Bbp+Aap))/2) - (((Ccp + Bbp))/2)*(((Ccp + Bbp))/2))*idz
        + gama*0.25*idz*(fabs(Bbp+Aap)*(Bbp-Aap)  -  fabs(Ccp+Bbp)*(Ccp-Bbp));

        return res;
}
