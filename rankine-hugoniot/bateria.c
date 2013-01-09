#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define METHODS 11
#define TESTS 2

int main()
{
	int i, j;
	char buffer[256];
	double U, T, A1[TESTS], A2[TESTS];
	const double B0[TESTS] =
	{
		0.,
		10.
	};
	const double Z[TESTS] =
	{
		2.,
		0.
	};
	const double S0[TESTS] =
	{
		0.02,
		0.1
	};
	const double r[TESTS] =
	{
		0.01,
		0.02
	};
	const double H1[TESTS] =
	{
		0.3,
		1.
	};
	const double Q1[TESTS] =
	{
		0.6671318441693031867,
		140.0175605713610490
	};
	const double H2[TESTS] =
	{
		0.2580357474523892850,
		0.03059155320377832031
	};
	const double Q2[TESTS] =
	{
		0.4463765053484111472,
		0.4712013382183945280
	};
	const double cfl[METHODS] =
	{
		0.9,
		0.9,
		4.,
		0.9,
		0.9,
		0.9,
		1.,
		0.75,
		4.,
		0.9,
		4.
	};
	const char *orden1[METHODS] =
	{
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 1 1 > i-1-1-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 2 1 > i-2-1-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 3 1 > i-3-1-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 4 1 > i-4-1-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 1 2 > i-1-2-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 2 2 > i-2-2-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 3 2 > i-3-2-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 1 3 > i-1-3-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 3 3 > i-3-3-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 1 4 > i-1-4-%d",
		"perl gen_bateria.pl %lg %lg %lg %lg %lg %lg %lg %lg %lg 3 4 > i-3-4-%d"
	};
	const char *orden2[METHODS] =
	{
		"../swocs i-1-1-%d o-1-1-%d f-1-1-%d a-1-1-%d\n",
		"../swocs i-2-1-%d o-2-1-%d f-2-1-%d a-2-1-%d\n",
		"../swocs i-3-1-%d o-3-1-%d f-3-1-%d a-3-1-%d\n",
		"../swocs i-4-1-%d o-4-1-%d f-4-1-%d a-4-1-%d\n",
		"../swocs i-1-2-%d o-1-2-%d f-1-2-%d a-1-2-%d\n",
		"../swocs i-2-2-%d o-2-2-%d f-2-2-%d a-2-2-%d\n",
		"../swocs i-3-2-%d o-3-2-%d f-3-2-%d a-3-2-%d\n",
		"../swocs i-1-3-%d o-1-3-%d f-1-3-%d a-1-3-%d\n",
		"../swocs i-3-3-%d o-3-3-%d f-3-3-%d a-3-3-%d\n",
		"../swocs i-1-4-%d o-1-4-%d f-1-4-%d a-1-4-%d\n",
		"../swocs i-3-4-%d o-3-4-%d f-3-4-%d a-3-4-%d\n"
	};
	const char *filename[METHODS + TESTS] =
	{
		"script-1-1",
		"script-2-1",
		"script-3-1",
		"script-4-1",
		"script-1-2",
		"script-2-2",
		"script-3-2",
		"script-1-3",
		"script-3-3",
		"script-1-4",
		"script-3-4",
		"s-1",
		"s-2"
	};
	FILE *file[METHODS + 1];
	for (j = 0; j < METHODS + TESTS; ++j) file[j] = fopen(filename[j], "w");
	for (i = 0; i < TESTS; ++i)
	{
		for (j = 0; j < METHODS; ++j)
		{
			snprintf(buffer, 256, orden1[j], B0[i], Z[i], S0[i], r[i],
				H1[i], Q1[i], H2[i], Q2[i], cfl[j], i + 1);
			system(buffer);
			fprintf(file[j], orden2[j], i + 1, i + 1, i + 1, i + 1);
		}
		A1[i] = H1[i] * (B0[i] + Z[i] * H1[i]);
		A2[i] = H2[i] * (B0[i] + Z[i] * H2[i]);
		U = (Q1[i] - Q2[i]) / (A1[i] - A2[i]);
		T = 50. / U;
		fprintf(file[METHODS + i], "0 %lg %lg 1 %lg %lg 1\n",
			A1[i], Q1[i], A1[i], Q1[i]);
		fprintf(file[METHODS + i], "25 %lg %lg 1 %lg %lg 1\n",
			A1[i], Q1[i], A1[i], Q1[i]);
		fprintf(file[METHODS + i], "25 %lg %lg 1 %lg %lg 0\n",
			A1[i], Q1[i], A2[i], Q2[i]);
		fprintf(file[METHODS + i], "%lg %lg %lg 1 %lg %lg 0\n",
			T * Q1[i] / A1[i], A1[i], Q1[i], A2[i], Q2[i]);
		fprintf(file[METHODS + i], "%lg %lg %lg 0 %lg %lg 0\n",
			T * Q1[i] / A1[i], A1[i], Q1[i], A2[i], Q2[i]);
		fprintf(file[METHODS + i], "75 %lg %lg 0 %lg %lg 1\n",
			A1[i], Q1[i], A2[i], Q2[i]);
		fprintf(file[METHODS + i], "75 %lg %lg 0 %lg %lg 0\n",
			A2[i], Q2[i], A2[i], Q2[i]);
		fprintf(file[METHODS + i], "100 %lg %lg 0 %lg %lg 0\n",
			A2[i], Q2[i], A2[i], Q2[i]);
	}
	for (j = 0; j < METHODS + TESTS; ++j) fclose(file[j]);
	system("chmod u+x script*");
	return 0;
}
