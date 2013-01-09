#define _GNU_SOURCE
#include <stdio.h>

int main(int argn, char **argc)
{
	int i, j, it, it2;
	double dr[3];
	char buffer[256];
	const char *label1[5] = {"1e-5", "1e-4", "1e-3", "1e-2", "1e-1"};
	const char *label2[5] = {"01", "02", "04", "08", "16"};
	FILE *r, *t, *t2, *s1[5], *s2[5];
	sprintf(buffer, "r-%s", argc[1]);
	r = fopen(buffer, "r");
	sprintf(buffer, "t-%s", argc[1]);
	t = fopen(buffer, "r");
	t2 = fopen("t-2-1", "r");
	for (i = 0; i < 5; ++i)
	{
		sprintf(buffer, "r-%s-%s", argc[1], label1[i]);
		s1[i] = fopen(buffer, "w");
		sprintf(buffer, "r-%s-%s", argc[1], label2[i]);
		s2[i] = fopen(buffer, "w");
	}
	for (i = 0; i < 5; ++i)
	{
		for (j = 0; j < 5; ++j)
		{
			fscanf(r, "%lf%lf%lf", dr, dr+1, dr+2);
			fscanf(t, "%d", &it);
			fscanf(t2, "%d", &it2);
			fprintf(s1[i], "%.6lg %.6lg %.6lg %d %d\n",
				dr[0], dr[1], dr[2], it, it2);
			fprintf(s2[j], "%.6lg %.6lg %.6lg %d %d\n",
				dr[0], dr[1], dr[2], it, it2);
		}
	}
	fcloseall();
	return 0;
}

