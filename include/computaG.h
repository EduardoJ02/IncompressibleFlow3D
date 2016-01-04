#ifndef COMPUTAG_INCLUDED
#define COMPUTAG_INCLUDED

#define B u3D[i][j][k]
#define C u3D[i-1][j][k]
#define D u3D[i][j+1][k]
#define Q u3D[i-1][j+1][k]

#define a v3D[i+1][j][k]
#define b v3D[i][j][k]
#define c v3D[i-1][j][k]
#define d v3D[i][j+1][k]
#define e v3D[i][j-1][k]

#define bfrente v3D[i][j][k+1]
#define batras v3D[i][j][k-1]
#define Bb w3D[i][j][k]
#define Aaj w3D[i][j+1][k]
#define Cc w3D[i][j][k-1]
#define Qqj w3D[i][j+1][k-1]

#define Esquerda 1
#define Superior 2
#define Inferior 4
#define Direita 8
#define Atras 16
#define Frente 32

#define Esquerdas flag3D[i][j][k]==Esquerda || flag3D[i][j][k]==Esquerda+Superior || flag3D[i][j][k]==Esquerda+Inferior || flag3D[i][j][k]==Esquerda+Atras || flag3D[i][j][k]==Esquerda+Superior+Atras || flag3D[i][j][k]==Esquerda+Inferior+Atras || flag3D[i][j][k]==Esquerda+Frente || flag3D[i][j][k]==Esquerda+Superior+Frente || flag3D[i][j][k]==Esquerda+Inferior+Frente
#define Direitas flag3D[i][j][k]==Direita || flag3D[i][j][k]==Direita+Superior || flag3D[i][j][k]==Direita+Inferior || flag3D[i][j][k]==Direita+Atras || flag3D[i][j][k]==Direita+Superior+Atras || flag3D[i][j][k]==Direita+Inferior+Atras || flag3D[i][j][k]==Direita+Frente || flag3D[i][j][k]==Direita+Superior+Frente || flag3D[i][j][k]==Direita+Inferior+Frente

#define Superiores flag3D[i][j][k]==Superior || flag3D[i][j][k]==Superior+Direita || flag3D[i][j][k]==Superior+Esquerda || flag3D[i][j][k]==Superior+Atras || flag3D[i][j][k]==Superior+Direita+Atras || flag3D[i][j][k]==Superior+Esquerda+Atras || flag3D[i][j][k]==Superior +Frente || flag3D[i][j][k]==Superior+Direita+Frente || flag3D[i][j][k]==Superior+Esquerda+Frente
#define Inferiores flag3D[i][j][k]==Inferior || flag3D[i][j][k]==Inferior+Direita || flag3D[i][j][k]==Inferior+Esquerda || flag3D[i][j][k]==Inferior+Atras || flag3D[i][j][k]==Inferior+Direita+Atras || flag3D[i][j][k]==Inferior+Esquerda+Atras || flag3D[i][j][k]==Inferior +Frente || flag3D[i][j][k]==Inferior+Direita+Frente || flag3D[i][j][k]==Inferior+Esquerda+Frente

#define Atrases    flag3D[i][j][k]==Atras || flag3D[i][j][k]==Atras+Inferior || flag3D[i][j][k]==Atras+Inferior+Direita || flag3D[i][j][k]==Atras+Inferior+Esquerda || flag3D[i][j][k]==Atras+Superior || flag3D[i][j][k]==Atras+Superior+Direita || flag3D[i][j][k]==Atras+Superior+Esquerda
#define Frentes    flag3D[i][j][k]==Frente || flag3D[i][j][k]==Frente+Inferior || flag3D[i][j][k]==Frente+Inferior+Direita || flag3D[i][j][k]==Frente+Inferior+Esquerda || flag3D[i][j][k]==Frente+Superior || flag3D[i][j][k]==Frente+Superior+Direita || flag3D[i][j][k]==Frente+Superior+Esquerda

extern double iRe, gama, dt, idx2, idx, idy2, idy, idz, idz2,  U_s;
extern double ***u3D, ***v3D;
extern int **flag;
extern int ***flag3D;

//double G(int i, int j, double **u, double **v);
//double Gc(int i, int j, double **u, double **v);
double G3D(int i, int j, int k, double ***u3D, double ***v3D);
double G3Dc(int i, int j, int k, double ***u3D, double ***v3D);
double G3DcExtra(int i, int j, int k, double ***u3D, double ***v3D,double ***w3D);

#endif
