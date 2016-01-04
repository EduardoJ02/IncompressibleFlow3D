#ifndef ARQUIVOS_H_INCLUDED
#define ARQUIVOS_H_INCLUDED

extern int iM; extern int jM; extern int ***flag3D; extern int **flag; extern int kM;

void Marcar(int **flag, char**flagchar,int iM, int jM);
char** leArquivo(void);
char** leArquivoDim(void);
char** escreveArquivo(double t, double ***u3D, double ***v3D, double ***w3D);

#endif
