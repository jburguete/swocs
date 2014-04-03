#include <stdio.h>
#include <math.h>

#define g 9.81L
#define n 0.02L
#define l 200.L
#define x0 50.L
#define xf 150.L
#define F "%.14Le"

long double b0, z;

long double a(long double h)
{return h * (b0 + z * h);}

long double u(long double s0, long double h)
{return sqrtl(s0) / n * powl(h, 2.L/3.L) * (b0 + 0.75L * z * h) / (b0 + z * h);}

long double q(long double s0, long double h)
{return a(h) * u(s0, h);}

long double v(long double s0, long double h1, long double h2)
{return (q(s0, h1) - q(s0, h2)) / (a(h1) - a(h2));}

long double beta(long double h)
{
	long double k = b0 + 0.75L * z * h;
	return 49.L/48.L * (b0 + z * h) * (b0 + 0.6L * z * h) / (k * k);
}

long double s0(long double s, long double h)
{return s * g * n * n / (beta(h) * powl(h, 1.L/3.L));}

long double t(long double s0, long double h1, long double h2)
{return (xf - x0) / v(s0, h1, h2);}

long double f1(long double x)
{
	long double k = (x - 1.L) / (powl(x, 2.L/3.L) - 1.L);
	return (x + 1.L) / (x + x) * k * k;
}

long double f2(long double x)
{
	long double k = x * (powl(x, 2.L/3.L) - 1.L);
	return 16.L/27.L * (x * x * x - 1.L) * (x * x - 1.L) / (k * k);
}

void write(char *filename, long double hl, long double hr,
	long double (*s)(long double), int scheme, int model)
{
	long double S, S0, Ql, Qr, Al, Ar;
	FILE *file;
	S = s(hl / hr);
	S0 = s0(S, hr);
	Ql = q(S0, hl);
	Qr = q(S0, hr);
	Al = a(hl);
	Ar = a(hr);
	file = fopen(filename, "w");
	fprintf(file, "2 2 2 1 1\n");
	fprintf(file, "2\n");
	fprintf(file, "0 "F" "F" "F" "F"\n", S0 * l, b0, z, 10.L + S0 * l);
	fprintf(file, F" 0 "F" "F" 10\n", l, b0, z);
	fprintf(file, F"\n", n);
	fprintf(file, "0 1 0 1\n");
	fprintf(file, "0\n");
	fprintf(file, "1\n");
	fprintf(file, "0 "F"\n", Ql);
	fprintf(file, "1\n");
	fprintf(file, "0 "F"\n", Ql);
	fprintf(file, "201 2\n");
	fprintf(file, "4\n");
	fprintf(file, "0 "F" "F" 1\n", Al, Ql);
	fprintf(file, F" "F" "F" 1\n", x0, Al, Ql);
	fprintf(file, F" "F" "F" 0\n", x0, Ar, Qr);
	fprintf(file, F" "F" "F" 0\n", l, Ar, Qr);
	fprintf(file, F" 0 0.9 0.01 %d 2 %d\n", t(S0, hl, hr), scheme, model);
	fclose(file);
}

int main()
{
	int i;
	int scheme[]={1, 4, 1, 1};
	int model[]={1, 1, 3, 4};
	char name[32];
	for (i = 0; i < 4; ++i)
	{
		b0 = 10.L;
		z = 0.L;
		snprintf(name, 32, "case1-%d-%d", scheme[i], model[i]);
		write(name, 1.1L, 1.L, &f1, scheme[i], model[i]);
		snprintf(name, 32, "case2-%d-%d", scheme[i], model[i]);
		write(name, 1L, 0.1L, &f1, scheme[i], model[i]);
		b0 = 0.L;
		z = 10.L;
		snprintf(name, 32, "case3-%d-%d", scheme[i], model[i]);
		write(name, 1.1L, 1.L, &f2, scheme[i], model[i]);
		snprintf(name, 32, "case4-%d-%d", scheme[i], model[i]);
		write(name, 1L, 0.2L, &f2, scheme[i], model[i]);
	}
	return 0;
}
