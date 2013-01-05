#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 200.
#define n 0.03
#define B 1.
#define METHODS 6

int main()
{
	int i, j;
	double T, Q, S, H;
	char buffer[256];
	const double t[25]=
	{
		1404.28, 743.629, 401.876, 220.482, 122.619,
		729.425, 397.856, 219.486, 122.381, 69.3709,
		358.492, 206.275, 118.156, 68.1019, 39.8050,
		143.001, 91.1419, 57.6769, 36.1538, 22.6192,
		47.5422, 32.1590, 22.7474, 16.4408, 11.7580
	};
	const double cfl[METHODS] =
	{
		0.9,
		0.9,
		4.,
		0.9,
		0.75,
		0.9
	};
	const char *orden1[METHODS] =
	{
		"perl gen_bateria.pl %lg %lg %lg %lg %lg 1 1 > i-1-1-%lg-%lg-%lg",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg 2 1 > i-2-1-%lg-%lg-%lg",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg 3 1 > i-3-1-%lg-%lg-%lg",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg 2 2 > i-2-2-%lg-%lg-%lg",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg 1 3 > i-1-3-%lg-%lg-%lg",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg 1 4 > i-1-4-%lg-%lg-%lg"
	};
	const char *orden2[METHODS] =
	{
		"../swocs i-1-1-%lg-%lg-%lg o-1-1-%lg-%lg-%lg "
		"f-1-1-%lg-%lg-%lg a-1-1-%lg-%lg-%lg\n",
		"../swocs i-2-1-%lg-%lg-%lg o-2-1-%lg-%lg-%lg "
		"f-2-1-%lg-%lg-%lg a-2-1-%lg-%lg-%lg\n",
		"../swocs i-3-1-%lg-%lg-%lg o-3-1-%lg-%lg-%lg "
		"f-3-1-%lg-%lg-%lg a-3-1-%lg-%lg-%lg\n",
		"../swocs i-2-2-%lg-%lg-%lg o-2-2-%lg-%lg-%lg "
		"f-2-2-%lg-%lg-%lg a-2-2-%lg-%lg-%lg\n",
		"../swocs i-1-3-%lg-%lg-%lg o-1-3-%lg-%lg-%lg "
		"f-1-3-%lg-%lg-%lg a-1-3-%lg-%lg-%lg\n",
		"../swocs i-1-4-%lg-%lg-%lg o-1-4-%lg-%lg-%lg "
		"f-1-4-%lg-%lg-%lg a-1-4-%lg-%lg-%lg\n"
	};
	const char *orden3[METHODS - 1] =
	{
		"./retraso %lg %lg %lg 2 1\n",
		"./retraso %lg %lg %lg 3 1\n",
		"./retraso %lg %lg %lg 2 2\n",
		"./retraso %lg %lg %lg 1 3\n",
		"./retraso %lg %lg %lg 1 4\n"
	};
	const char *filename[2 * METHODS - 1] =
	{
		"script-1-1",
		"script-2-1",
		"script-3-1",
		"script-2-2",
		"script-1-3",
		"script-1-4",
		"script-r-2-1",
		"script-r-3-1",
		"script-r-2-2",
		"script-r-1-3",
		"script-r-1-4"
	};
	FILE *file[2 * METHODS - 1];
	for (j = 0; j < 2 * METHODS - 1; ++j) file[j] = fopen(filename[j], "w");
	for (S = 0.00001, i = 0; S < 0.2; S *= 10.)
	{
		for (H = 0.1; H < 2.; H *= 2., ++i)
		{
			Q = sqrt(S) * pow(H, 5./3.) / n;
			T = t[i];
			for (j = 0; j < METHODS; ++j)
			{
				snprintf(buffer, 256, orden1[j], L, S, Q, T, cfl[j], S, H, T);
				system(buffer);
				fprintf(file[j], orden2[j], S, H, T, S, H, T, S, H, T, S, H, T);
			}
			for (j = 0; j < METHODS - 1; ++j)
				fprintf(file[METHODS + j], orden3[j], S, H, T);
		}
	}
	for (j = 0; j < 2 * METHODS - 1; ++j) fclose(file[j]);
	system("chmod u+x script*");
	return 0;
}
