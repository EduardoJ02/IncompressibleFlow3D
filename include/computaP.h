#ifndef COMPUTAP_H_INCLUDED
#define COMPUTAP_H_INCLUDED

#include "Alocar.h"

#define PDirC mP_ant3D[i][j][k]-dx*idt*(u3D[i+1][j][k]- mF3D[i+1][j][k])
#define PEsqC mP_ant3D[i][j][k]+dx*idt*(u3D[i][j][k]- mF3D[i][j][k])
#define PInfC mP_ant3D[i][j][k]+dy*idt*(v3D[i][j][k]- mG3D[i][j][k])
#define PSupC mP_ant3D[i][j][k]-dy*idt*(v3D[i][j+1][k]- mG3D[i][j+1][k])

#define PAtrasC mP_ant3D[i][j][k]+dz*idt*(w3D[i][j][k]- mH3D[i][j][k])
#define PFrenteC mP_ant3D[i][j][k]-dz*idt*(w3D[i][j][k+1]- mH3D[i][j][k+1])

#define PSup mP3D[i][j+1][k]
#define PInf mP_ant3D[i][j-1][k]
#define PDir mP3D[i+1][j][k]
#define PEsq mP_ant3D[i-1][j][k]

#define PFrente mP3D[i][j][k+1]
#define PAtras mP3D[i][j][k-1]

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


extern int ***flag3D;
extern int iM, jM, kM;
extern double idx2, idy2, dt, dx, dy, dz, idz2, idt;

void computaP(double ***mP3D, double ***mdP3D, double ***mP_ant3D, double ***residual3D, double ***mF3D, double ***mG3D, double ***mH3D, double ***u3D, double ***v3D, double ***w3D, int k);

#endif
