#include <stdio.h>

int main(int argn, char **argc)
{
	FILE *input_file, *output_file;
	double x[5];
	int i[4];
	if (argn != 3)
	{
		printf("The sintaxis is:\ntranslate-1-2 input_file output_file\n");
		return 1;
	}
	input_file = fopen(argc[1], "r");
	if (!input_file)
	{
		printf("Unable to open the input file\n");
		return 2;
	}
	if (fscanf(input_file, "%lf%lf%lf%lf%lf%d%d%d%d\n",
		x, x+1, x+2, x+3, x+4, i, i+1, i+2, i+3) != 9)
	{
		printf("Bad input file\n");
		return 3;
	}
	fprintf(output_file, "%lg %lg %lg %lg %d %d %d %d\n" 

}
