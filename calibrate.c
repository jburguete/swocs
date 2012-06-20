#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PESO_AVANCE 0.5
#define PESO_CALADO 0.1
#define PESO_CONCENTRACION 0.1

double file_mean_square_error
	(char *namea,int ixa,int ifa,int na,char *namer,int ixr,int ifr,int nr)
{
	int i,j,endr;
	double k,k2,xa,fa,xr1,fr1,xr2,fr2,aa[na],ar[nr];
	FILE *filea,*filer;
	endr=i=0;
	k=0.;
	filea=fopen(namea,"r");
	if (!filea) goto exit_mse;
	filer=fopen(namer,"r");
	if (!filer) goto exit_mse;
	for (j=0; j<nr; ++j)
		if (fscanf(filer,"%lf",ar+j)!=1) goto exit_mse;
	xr1=ar[ixr-1];
	fr1=ar[ifr-1];
	for (j=0; j<nr; ++j)
		if (fscanf(filer,"%lf",ar+j)!=1) endr=1;
	xr2=ar[ixr-1];
	fr2=ar[ifr-1];
	for (i=0; !endr; ++i)
	{
		for (j=0; j<na; ++j)
			if (fscanf(filea,"%lf",aa+j)!=1) goto exit_mse;
		xa=aa[ixa-1];
		fa=aa[ifa-1];
		while (xa>xr2)
		{
			xr1=xr2;
			fr1=fr2;
			for (j=0; j<nr; ++j) if (fscanf(filer,"%lf",ar+j)!=1)
			{
				endr=1;
				goto end_filer;
			}
			xr2=ar[ixr-1];
			fr2=ar[ifr-1];
		}
end_filer:
		if (!endr && xa>xr1)
			k2=fa-fr1-(xa-xr1)*(fr2-fr1)/(xr2-xr1);
		else k2=fa-fr1;
		k+=k2*k2;
	}
	for (; 1; ++i)
	{
		for (j=0; j<na; ++j)
			if (fscanf(filea,"%lf",aa+j)!=1) goto exit_mse;
		xa=aa[ixa-1];
		fa=aa[ifa-1];
		k2=fa-fr1;
		k+=k2*k2;
	}
exit_mse:
	if (i==0) return 0.;
	return sqrt(k/i);
}

int main(int argn, char **argc)
{
	int i,j;
	double etotal,e[4][6];

	// simulamos los 4 casos
	system("./script_nery");

	// errores en avance
	e[0][0]=file_mean_square_error("nery/q1a",1,2,2,"nery/advances1.out",2,1,2);
	e[1][0]=file_mean_square_error("nery/q2a",1,2,2,"nery/advances2.out",2,1,2);
	e[2][0]=file_mean_square_error("nery/q3a",1,2,2,"nery/advances3.out",2,1,2);
	e[3][0]=file_mean_square_error("nery/q4a",1,2,2,"nery/advances4.out",2,1,2);

	// errores en calados aguas arriba
	e[0][1]=file_mean_square_error("nery/h1",1,2,2,"nery/probes1.out",1,2,11);
	e[1][1]=file_mean_square_error("nery/h2",1,2,2,"nery/probes2.out",1,2,11);
	e[2][1]=file_mean_square_error("nery/h3",1,2,2,"nery/probes3.out",1,2,11);
	e[3][1]=file_mean_square_error("nery/h4",1,2,2,"nery/probes4.out",1,2,11);

	// errores en sondas de concentración
	e[0][2]=file_mean_square_error("nery/q1c1",1,3,3,"nery/probes1.out",1,5,11);
	e[1][2]=file_mean_square_error("nery/q2c1",1,3,3,"nery/probes2.out",1,5,11);
	e[2][2]=file_mean_square_error("nery/q3c1",1,3,3,"nery/probes3.out",1,5,11);
	e[3][2]=file_mean_square_error("nery/q4c1",1,3,3,"nery/probes4.out",1,5,11);
	e[0][3]=file_mean_square_error("nery/q1c2",1,3,3,"nery/probes1.out",1,7,11);
	e[1][3]=file_mean_square_error("nery/q2c2",1,3,3,"nery/probes2.out",1,7,11);
	e[2][3]=file_mean_square_error("nery/q3c2",1,3,3,"nery/probes3.out",1,7,11);
	e[3][3]=file_mean_square_error("nery/q4c2",1,3,3,"nery/probes4.out",1,7,11);
	e[0][4]=file_mean_square_error("nery/q1c3",1,3,3,"nery/probes1.out",1,9,11);
	e[1][4]=file_mean_square_error("nery/q2c3",1,3,3,"nery/probes2.out",1,9,11);
	e[2][4]=file_mean_square_error("nery/q3c3",1,3,3,"nery/probes3.out",1,9,11);
	e[3][4]=file_mean_square_error("nery/q4c3",1,3,3,"nery/probes4.out",1,9,11);
	e[0][5]=file_mean_square_error("nery/q1c4",1,3,3,"nery/probes1.out",1,11,11);
	e[1][5]=file_mean_square_error("nery/q2c4",1,3,3,"nery/probes2.out",1,11,11);
	e[2][5]=file_mean_square_error("nery/q3c4",1,3,3,"nery/probes3.out",1,11,11);
	e[3][5]=file_mean_square_error("nery/q4c4",1,3,3,"nery/probes4.out",1,11,11);

	// normalizamos avances
	double a[4]={2747.,1874.,1280.,978.};
	for (i=0; i<4; ++i) e[i][0]/=a[i];

	// normalizamos calados
	double h[4]={0.0744,0.0703,0.111,0.0998};
	for (i=0; i<4; ++i) e[i][1]/=h[i];

	// normalizamos concentraciones
	for (j=2; j<6; ++j)
		for (i=0; i<4; ++i) e[i][j]/=10.6;

/*
for (j=0; j<6; ++j)
	for (i=0; i<4; ++i)
		printf("i=%d j=%d e=%lg\n",i,j,e[i][j]);
*/

	// pesamos los errores
	etotal=0.;
	for (i=0; i<4; ++i) etotal+=PESO_AVANCE*e[i][0]*e[i][0];
	for (i=0; i<4; ++i) etotal+=PESO_CALADO*e[i][1]*e[i][1];
	for (j=2; j<6; ++j)
		for (i=0; i<4; ++i) etotal+=PESO_CONCENTRACION*e[i][j]*e[i][j];

	// escribimos en pantalla
	printf("%lg\n",sqrt(etotal/4.));
	return 0;
}
