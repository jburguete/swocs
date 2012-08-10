#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 200.
#define n 0.03
#define B 1.

int main()
{
	double T, Q, S;
	char buffer[256];
	FILE *file1, *file2, *file3, *file4, *file5, *file6, *file7;
	file1 = fopen("script-2-1", "w");
	file2 = fopen("script-2-2", "w");
	file3 = fopen("script-1-3", "w");
	file4 = fopen("script-1-4", "w");
	file5 = fopen("script-r-2-2", "w");
	file6 = fopen("script-r-1-3", "w");
	file7 = fopen("script-r-1-4", "w");
	for (S = 0.00001; S < 0.2; S *= 10.)
	{
		for (Q = 0.001; Q < 20.; Q *= 10.)
		{
			T = round(2.1 * L * pow(B / Q, 0.4) * pow(n, 0.6) * pow(S, -0.1));
			snprintf(buffer, 256,
				"perl gen_bateria.pl %lg %lg %lg %lg 2 1 > i-2-1-%lg-%lg-%lg",
				L, S, Q, T, S, Q, T);
			system(buffer);
			fprintf(file1, "../swocs i-2-1-%lg-%lg-%lg o-2-1-%lg-%lg-%lg "
				"f-2-1-%lg-%lg-%lg a-2-1-%lg-%lg-%lg\n",
				S, Q, T, S, Q, T, S, Q, T, S, Q, T);
			snprintf(buffer, 256,
				"perl gen_bateria.pl %lg %lg %lg %lg 2 2 > i-2-2-%lg-%lg-%lg",
				L, S, Q, T, S, Q, T);
			system(buffer);
			fprintf(file2, "../swocs i-2-2-%lg-%lg-%lg o-2-2-%lg-%lg-%lg "
				"f-2-2-%lg-%lg-%lg a-2-2-%lg-%lg-%lg\n",
				S, Q, T, S, Q, T, S, Q, T, S, Q, T);
			snprintf(buffer, 256,
				"perl gen_bateria.pl %lg %lg %lg %lg 1 3 > i-1-3-%lg-%lg-%lg",
				L, S, Q, T, S, Q, T);
			system(buffer);
			fprintf(file3, "../swocs i-1-3-%lg-%lg-%lg o-1-3-%lg-%lg-%lg "
				"f-1-3-%lg-%lg-%lg a-1-3-%lg-%lg-%lg\n",
				S, Q, T, S, Q, T, S, Q, T, S, Q, T);
			snprintf(buffer, 256,
				"perl gen_bateria.pl %lg %lg %lg %lg 1 4 > i-1-4-%lg-%lg-%lg",
				L, S, Q, T, S, Q, T);
			system(buffer);
			fprintf(file4, "../swocs i-1-4-%lg-%lg-%lg o-1-4-%lg-%lg-%lg "
				"f-1-4-%lg-%lg-%lg a-1-4-%lg-%lg-%lg\n",
				S, Q, T, S, Q, T, S, Q, T, S, Q, T);
			fprintf(file5, "./retraso %lg %lg %lg 2 1 2 2\n", S, Q, T);
			fprintf(file6, "./retraso %lg %lg %lg 2 1 1 3\n", S, Q, T);
			fprintf(file7, "./retraso %lg %lg %lg 2 1 1 4\n", S, Q, T);
		}
	}
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	fclose(file5);
	fclose(file6);
	fclose(file7);
	return 0;
}
