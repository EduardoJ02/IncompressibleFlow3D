#ifndef DISCRET_H_INCLUDED
#define DISCRET_H_INCLUDED

    double u2_x2(int i, int j, double Dp, double Bp, double Ep);

    double u2_y2(int i, int j, double Ap, double Bp, double Cp);

    double u2(int i, int j, double Bp, double Dp, double Ep);

    double uv_y(int i, int j, double bp, double dp, double Bp, double Ap, double cp, double qp, double Cp);

    //Extra
    double u2_z2(int i, int j, double Bfrentep, double Bp, double Batrasp);
    double uw_z(int i, int j, double Bbp, double Aaip,double Bp,double Bfrentep,double Ccp,double Qqip,double Batrasp);


    double v2_x2(int i, int j, double dp, double bp, double ep);

    double v2_y2(int i, int j, double ap, double bp, double cp);

    double v2(int i, int j, double bp, double ap, double cp);

    double uv_x(int i, int j, double Bp, double Ap, double bp, double dp, double Ep, double Qp, double ep);

    //Extra
    double v2_z2(int i, int j, double bfrentep, double bp, double batrasp);
    double vw_z(int i, int j, double Bbp, double Aajp,double bp,double bfrentep,double Ccp,double Qqjp,double batrasp);



    double laplaceW(int i, int j, double Aaip,double Bbp,double Ccip,double Aajp,double Ccjp,double Aap,double Ccp);
    double uw_x(int i, int j, double Bp,double Bfp,double Bbp,double Aaip,double Cp,double Cfp,double Ccip);
    double vw_y(int i, int j, double bp,double bfp,double Bbp,double Aajp,double ep,double efp,double Ccjp);
    double w2(int i, int j, double Bbp,double Aap,double Ccp);

#endif
