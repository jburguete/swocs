#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#define g 9.81L
#define n 0.02L
#define l 1.L
#define lambda 0.1L
#define k 2.L * M_PIl / lambda
#define F "%.14Le"
#define beta 49.L/48.L
#define N2 2001
#define cfl 0.9L

long double b0, z;

long double a(long double h)
{return h * (b0 + z * h);}

long double b(long double h)
{return b0 + 2.L * z * h;}

long double c02(long double h0)
{return g * a(h0) / b(h0);}

long double c0(long double h0)
{return sqrtl(c02(h0));}

long double s0(long double h0, long double u0)
{
	long double nu;
	nu = n * u0 * (b0 + z * h0) / (b0 + 0.75L * z * h0);
	return nu * nu * pow(h0, -4.L/3.L);
}

long double v1(long double h0, long double u0)
{
	return beta * u0 + sqrtl(c02(h0) + beta * (beta - 1.L) * u0 * u0);
}

long double v2(long double h0, long double u0)
{
	return beta * u0 - sqrtl(c02(h0) + beta * (beta - 1.L) * u0 * u0);
}

long double u1(long double h0, long double dh, long double u0, long double v)
{
	return (v - u0) * dh / h0;
}

long double tf(long double h0, long double u0)
{return 0.5L * lambda / v1(h0, u0);}

long double dh(long double h1, long double v, long double phi, long double x,
	long double t)
{
	return h1 * sinl(k * (x - v * t) + phi);
}

long double h(long double h0, long double h1, long double h2, long double u0,
	long double phi, long double x, long double t)
{
	return h0 + dh(h1, v1(h0, u0), phi, x, t)
		+ dh(h2, v2(h0, u0), phi, x, t);
}

long double du(long double u1, long double v, long double phi, long double x,
	long double t)
{
	return u1 * sinl(k * (x - v * t) + phi);
}

long double u(long double h0, long double h1, long double h2, long double u0,
	long double phi, long double x, long double t)
{
	long double U1, U2, V1, V2;
	V1 = v1(h0, u0);
	V2 = v2(h0, u0);
	U1 = u1(h0, h1, u0, V1);
	U2 = u1(h0, h2, u0, V2);
	return u0 + du(U1, V1, phi, x, t) + du(U2, V2, phi, x, t);
}

long double q(long double h0, long double h1, long double h2, long double u0,
	long double phi, long double x, long double t)
{return a(h(h0, h1, h2, u0, phi, x, t)) * u(h0, h1, h2, u0, phi, x, t);}

void write(char *filename1, char *filename2, long double h0, long double h1,
	long double h2, long double u0, long double phi, int N, int scheme,
	int model)
{
	int i, NT;
	long double S0, TF, X, T;
	FILE *file;
	S0 = s0(h0, u0);
	TF = tf(h0, u0);
	printf("v1=%Lg v2=%Lg tf=%Lg\n", v1(h0, u0), v2(h0, u0), TF);
	NT = ceill(0.5L * lambda * N / (cfl * l)); 
	file = fopen(filename1, "w");
	fprintf(file, "2 2 2 1 1\n");
	fprintf(file, "2\n");
	fprintf(file, "0 "F" "F" "F" "F"\n", S0 * l, b0, z, 10.L + S0 * l);
	fprintf(file, F" 0 "F" "F" 10\n", l, b0, z);
	fprintf(file, F"\n", n);
	fprintf(file, "0 1 0 1\n");
	fprintf(file, "0\n");
	fprintf(file, "%d\n", NT);
	for (i = 0; i < NT; ++i)
	{
		T = i * TF / (NT - 1);
		fprintf(file, "0 "F"\n", q(h0, h1, h2, u0, phi, 0.L, T));
	}
	fprintf(file, "1\n");
	fprintf(file, "0 0\n");
	fprintf(file, "%d 2\n", N);
	fprintf(file, "%d\n", N);
	for (i = 0; i < N; ++i)
	{
		X = i * l / (N - 1);
		fprintf(file, F" "F" "F" 0\n", X, a(h(h0, h1, h2, u0, phi, X, 0.L)),
			q(h0, h1, h2, u0, phi, X, 0.L));
	}
	fprintf(file, F" 0 0.9 0.01 %d 2 %d\n", TF, scheme, model);
	fclose(file);
	file = fopen(filename2, "w");
	for (i = 0; i < N2; ++i)
	{
		X = i * l / (N2 - 1);
		fprintf(file, F" "F"\n", X, h(h0, h1, h2, u0, phi, X, TF));
	}
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
		write(name, "sol1", 1.L, 1e-7L, 1e-7L, 1.L, 0.L, 201, scheme[i],
			model[i]);
		snprintf(name, 32, "case2-%d-%d", scheme[i], model[i]);
		write(name, "sol2", 1.L, 1e-7L, 1e-7L, 4.L, 0.L, 201, scheme[i],
			model[i]);
		b0 = 0.L;
		z = 10.L;
		snprintf(name, 32, "case3-%d-%d", scheme[i], model[i]);
		write(name, "sol3", 1.L, 1e-7L, 1e-7L, 1.L, 0.L, 201, scheme[i],
			model[i]);
		snprintf(name, 32, "case4-%d-%d", scheme[i], model[i]);
		write(name, "sol4", 1.L, 1e-7L, 1e-7L, 3.L, 0.L, 201, scheme[i],
			model[i]);
		snprintf(name, 32, "case4b-%d-%d", scheme[i], model[i]);
		write(name, "sol4b", 1.L, 1e-7L, 1e-7L, 3.L, 0.L, 801, scheme[i],
			model[i]);
	}
	return 0;
}
