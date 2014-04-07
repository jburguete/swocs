#include <stdio.h>
#include <math.h>

#define g 9.81L
#define l 200.L
#define x0 50.L
#define xf 150.L
#define F "%.14Le"

long double b0, z;

long double a(long double h)
{return h * (b0 + z * h);}

long double i1(long double h)
{return h * h * (0.5L * b0 + 1.L/3.L * z * h);}

long double beta(long double h)
{
	long double k = b0 + 0.75L * z * h;
	return 49.L/48.L * (b0 + z * h) * (b0 + 0.6L * z * h) / (k * k);
}

void write(char *filename, long double hl, long double hr, long double ur,
	int n, long double cfl, int outlet, int scheme, int model)
{
	long double ul, v, Ql, Qr, Al, Ar;
	FILE *file;
	Al = a(hl);
	Ar = a(hr);
	Qr = Ar * ur;
	ul = ur + sqrtl(g * (i1(hl) - i1(hr)) * (Al - Ar) / (beta(hr) * Al * Ar));
	Ql = Al * ul;
	v = (Ql - Qr) / (Al - Ar);
	file = fopen(filename, "w");
	fprintf(file, "2 %d 2 1 1\n", outlet);
	fprintf(file, "2\n");
	fprintf(file, "0 0 "F" "F" 10\n", b0, z);
	fprintf(file, F" 0 "F" "F" 10\n", l, b0, z);
	fprintf(file, "0\n");
	fprintf(file, "0 1 0 1\n");
	fprintf(file, "0\n");
	fprintf(file, "1\n");
	fprintf(file, "0 "F"\n", Ql);
	fprintf(file, "1\n");
	fprintf(file, "0 "F"\n", Ql);
	fprintf(file, "%d 2\n", n);
	fprintf(file, "4\n");
	fprintf(file, "0 "F" "F" 1\n", Al, Ql);
	fprintf(file, F" "F" "F" 1\n", x0, Al, Ql);
	fprintf(file, F" "F" "F" 0\n", x0, Ar, Qr);
	fprintf(file, F" "F" "F" 0\n", l, Ar, Qr);
	fprintf(file, F" 0 "F" 0.01 %d 2 %d\n", (xf - x0) / v, cfl, scheme, model);
	fclose(file);
}

int main()
{
	int i;
	int scheme[]={1, 4};
	int model[]={1, 1};
	char name[32];
	for (i = 0; i < 2; ++i)
	{
		b0 = 10.L;
		z = 0.L;
		snprintf(name, 32, "case1-%d-%d", scheme[i], model[i]);
		write(name, 1.1L, 1.L, 0.L, 201, 0.9, 1, scheme[i], model[i]);
		snprintf(name, 32, "case2-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 0.1L, 0.L, 201, 0.9, 1, scheme[i], model[i]);
		snprintf(name, 32, "case3-%d-%d", scheme[i], model[i]);
		write(name, 1.1L, 1.L, 4.L, 201, 0.9, 2, scheme[i], model[i]);
		snprintf(name, 32, "case4-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 0.1L, 1.L, 201, 0.9, 2, scheme[i], model[i]);
		b0 = 0.L;
		z = 10.L;
		snprintf(name, 32, "case5-%d-%d", scheme[i], model[i]);
		write(name, 1.1L, 1.L, 0.L, 201, 0.9, 1, scheme[i], model[i]);
		snprintf(name, 32, "case6-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 0.2L, 0.L, 201, 0.9, 1, scheme[i], model[i]);
		snprintf(name, 32, "case7-%d-%d", scheme[i], model[i]);
		write(name, 1.1L, 1.L, 3.L, 201, 0.9, 2, scheme[i], model[i]);
		snprintf(name, 32, "case8-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 0.2L, 1.L, 201, 0.9, 2, scheme[i], model[i]);
		snprintf(name, 32, "case8b-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 0.2L, 1.L, 401, 0.9, 2, scheme[i], model[i]);
		snprintf(name, 32, "case8c-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 0.2L, 1.L, 201, 0.1, 2, scheme[i], model[i]);
	}
	return 0;
}
