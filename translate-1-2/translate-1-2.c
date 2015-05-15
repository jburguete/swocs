#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argn, char **argc)
{
	FILE *input_file, *output_file;
	double x[5];
	int i[4];
	size_t n = 0;
	char *buffer = NULL;
	if (argn != 3)
	{
		printf("The syntax is:\ntranslate-1-2 input_file output_file\n");
		return 1;
	}
	input_file = fopen(argc[1], "r");
	if (!input_file)
	{
		printf("Unable to open the input file\n");
		return 2;
	}
	output_file = fopen(argc[2], "w");
	if (!output_file)
	{
		printf("Unable to open the output file\n");
		return 2;
	}
	if (fscanf(input_file, "%lf%lf%lf%lf%lf%d%d%d%d\n",
		x, x+1, x+2, x+3, x+4, i, i+1, i+2, i+3) != 9)
	{
		printf("Bad input file\n");
		return 3;
	}
	fprintf(output_file, "%lg %lg %lg %lg\n%d\n%d %d %d\n2\n0 %lg\n%lg 0\n",
		x[0], x[2], x[3], x[4], i[0], i[1], i[2], i[3], x[1] * x[0], x[0]);
	while (getline(&buffer, &n, input_file) >= 0)
	{
		fprintf(output_file, buffer);
		n = 0;
		if (!strcmp(buffer,
			"(friction model) (infiltration model) (diffusion model)\n"))
				fprintf(output_file, "(points number of geometry)\n"
					"(x-coordinate of the geometry point) "
					"(z-coordinate of the geometry point)\n"
					"...\n");
	}
	fclose(input_file);
	fclose(output_file);
	return 0;
}
