#include <stdio.h>

int main(int argn, char **argc)
{
	int a[2];
	double x[2], y[1], S, Q, T;
	char buffer[256];
	FILE *file;
	if (argn != 6) return 1;
	sscanf(argc[1], "%lf", &S);
	sscanf(argc[2], "%lf", &Q);
	sscanf(argc[3], "%lf", &T);
	sscanf(argc[4], "%u", a);
	sscanf(argc[5], "%u", a+1);
	snprintf(buffer, 256, "a-%u-%u-%lg-%lg-%lg", a[0], a[1], S, Q, T);
	file = fopen(buffer, "r");
	if (!file)
	{
		printf("Unable to open the file %s\n", buffer);
		return 2;
	}
	while (fscanf(file, "%lf%lf", x, x+1)==2) y[0]=x[1];
	printf("Delay=%lg\n", y[0]);
	fclose(file);
	snprintf(buffer, 256, "r-%u-%u", a[0], a[1]);
	file = fopen(buffer, "a");
	fprintf(file, "%lf %lf %lf\n", S, Q, y[0]);
	fclose(file);
	return 0;
}
