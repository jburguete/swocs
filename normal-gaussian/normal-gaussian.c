#include <stdio.h>
#include <math.h>

#define g 9.81L
#define n 0.02L
#define l 200.L
#define x0 50.L
#define xf 150.L
#define F "%.14Le"
#define beta 49.L/48.L
#define N 201
#define b0 10.L

long double a(long double h)
{return h * b0;}

long double c02(long double h0)
{return g * h0;}

long double c0(long double h0)
{return sqrtl(c02(h0));}

long double u0(long double h0)
{return 3.L * c0(h0) / sqrtl(27.L - 21.L * beta);}

long double s0(long double h0)
{
	long double nu;
	nu = n * u0(h0);
	return nu * nu * pow(h0, -4.L/3.L);
}

long double q0(long double h0)
{return a(h0) * u0(h0);}

long double v0(long double h0)
{
	long double U;
	U = u0(h0);
	return beta * U + sqrtl(c02(h0) + beta * (beta - 1.L) * U * U);
}

long double du(long double h0, long double dh)
{
	return (v0(h0) - u0(h0)) * dh / h0;
}

long double t(long double h0)
{return (xf - x0) / v0(h0);}

long double h(long double h0, long double dh, long double x, long double width)
{
	long double k;
	k = (x - x0) / width;
	return h0 + dh * expl(-k * k);
}

long double u(long double u0, long double du, long double x, long double width)
{
	long double k;
	k = (x - x0) / width;
	return u0 + du * expl(-k * k);
}

long double q(long double h0, long double dh, long double u0, long double du,
	long double x, long double width)
{
	return a(h(h0, dh, x, width)) * u(u0, du, x, width);
}

void write(char *filename, long double h0, long double dh, long double width,
	int scheme, int model)
{
	int i;
	long double S0, H0, U0, DU, x;
	FILE *file;
	S0 = s0(h0);
	H0 = h(h0, dh, 0, width);
	U0 = u0(h0);
	DU = du(h0, dh);
	file = fopen(filename, "w");
	fprintf(file, "2 2 2 1 1\n");
	fprintf(file, "2\n");
	fprintf(file, "0 "F" "F" 0 "F"\n", S0 * l, b0, 10.L + S0 * l);
	fprintf(file, F" 0 "F" 0 10\n", l, b0);
	fprintf(file, F"\n", n);
	fprintf(file, "0 1 0 1\n");
	fprintf(file, "0\n");
	fprintf(file, "1\n");
	fprintf(file, "0 "F"\n", q(H0, 0.L, U0, 0.L, 0.L, width));
	fprintf(file, "1\n");
	fprintf(file, "0 0\n");
	fprintf(file, "%d 2\n", N);
	fprintf(file, "%d\n", N);
	for (i = 0; i < N; ++i)
	{
		x = i * l / (N - 1);
		H0 = h(h0, dh, x, width);
		fprintf(file, F" "F" "F" 0\n", x, a(H0), q(h0, dh, U0, DU, x,width));
	}
	fprintf(file, F" 0 0.9 0.01 %d 2 %d\n", t(h0), scheme, model);
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
		snprintf(name, 32, "case1-%d-%d", scheme[i], model[i]);
		write(name, 1.L, 1e-7L, 10.L, scheme[i], model[i]);
	}
	return 0;
}
