#include <stdio.h>

#define L 100.

int main(int argn, char **argc)
{
	int a[4];
	double x[2], y[4], S, Q, T;
	char buffer[256];
	FILE *file;
	if (argn != 8) return 1;
	sscanf(argc[1], "%lf", &S);
	sscanf(argc[2], "%lf", &Q);
	sscanf(argc[3], "%lf", &T);
	sscanf(argc[4], "%u", a);
	sscanf(argc[5], "%u", a+1);
	sscanf(argc[6], "%u", a+2);
	sscanf(argc[7], "%u", a+3);
	snprintf(buffer, 256, "a-%u-%u-%lg-%lg-%lg", a[0], a[1], S, Q, T);
	file = fopen(buffer, "r");
	if (!file)
	{
		printf("Unable to open the file %s\n", buffer);
		return 2;
	}
	while (fscanf(file, "%lf%lf", x, x+1)==2) if (x[1]==100.) break;
	printf("x=%lg\n", x[1]);
	if (x[1]!=100.) return 2;
	fclose(file);
	snprintf(buffer, 256, "a-%u-%u-%lg-%lg-%lg", a[2], a[3], S, Q, T);
	file = fopen(buffer, "r");
	if (!file)
	{
		printf("Unable to open the file %s\n", buffer);
		return 2;
	}
	for (; fscanf(file, "%lf%lf", y, y+1)==2; y[2]=y[0], y[3]=y[1])
		if (y[0]>x[0]) break;
	fclose(file);
	printf("x %lf %lf\n", x[0], x[1]);
	printf("y %lf %lf %lf %lf\n", y[0], y[1], y[2], y[3]);
	x[0]=y[1]+(y[3]-y[1])*(x[0]-y[2])/(y[0]-y[2])-x[1];
	printf("Delay=%lf\n", x[0]);
	snprintf(buffer, 256, "r-%u-%u", a[2], a[3]);
	file = fopen(buffer, "a");
	fprintf(file, "%lf %lf %lf\n", S, Q, x[0]);
	fclose(file);
	return 0;
}
