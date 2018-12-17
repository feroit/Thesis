#include <stdlib.h>
#include <cmath>
#define M_PI (3.14159265358979323846)

typedef double function_callback(double phi, double psi);
function_callback* u0;
typedef double function_callback2(double phi, double psi, double t);
function_callback2* f;

int inline g(int i, int size)
{
	if (i < 0) i = size - 1;
	if (i > size) i = 1;
	return i;
}
/*double inline u0(double phi, double psi)
{
	return 8.0*cos(4 * psi)*cos(4 * phi);
}
double inline f(double phi, double psi, double t)
{
	return 0.0;
}*/
void inline reduce(double* mas, double* A, double* C, double* B, double* F, int S)
{
	double *p = new double[S], *q = new double[S];
	double *alpha = new double[S + 2], *beta = new double[S + 2], *gamma = new double[S + 2];
	alpha[2] = B[1] / C[1];
	beta[2] = F[1] / C[1];
	gamma[2] = A[1] / C[1];
	for (int i = 2; i <= S; i++)
	{
		alpha[i + 1] = B[i] / (C[i] - alpha[i] * A[i]);
		beta[i + 1] = (F[i] + beta[i] * A[i]) / (C[i] - alpha[i] * A[i]);
		gamma[i + 1] = (A[i] * gamma[i]) / (C[i] - alpha[i] * A[i]);
	}
	p[S - 1] = beta[S];
	q[S - 1] = alpha[S] + gamma[S];
	for (int i = S - 2; i>0; i--)
	{
		p[i] = alpha[i + 1] * p[i + 1] + beta[i + 1];
		q[i] = alpha[i + 1] * q[i + 1] + gamma[i + 1];
	}
	mas[S] = (p[1] * alpha[S + 1] + beta[S + 1]) / (1 - gamma[S + 1] - q[1] * alpha[S + 1]);
	mas[0] = mas[S];
	for (int i = 1; i<S; i++)
	{
		mas[i] = p[i] + mas[S] * q[i];
	}
	delete[] p;
	delete[] q;
	delete[] alpha;
	delete[] beta;
	delete[] gamma;
}
__declspec(dllexport) void inline Solution(double** data, function_callback fcb, function_callback2 fcb2, double* param, int* size)
{
	u0 = fcb;
	f = fcb2;
	int T = size[0];
	int N = size[1];
	int M = size[2];
	double ***u = new double **[T];
	double ***fmas = new double **[T];
	for (int k = 0; k<T; k++)
	{
		u[k] = new double *[N + 1];
		fmas[k] = new double *[N + 1];
		for (int i = 0; i < N + 1; i++)
		{
			u[k][i] = new double[M + 1];
			fmas[k][i] = new double[M + 1];
		}
	}
	double **u_half = new double *[N + 1];
	for (int j = 0; j<N + 1; j++)
		u_half[j] = new double[M + 1];
	double *AM = new double[M + 1], *CM = new double[M + 1], *BM = new double[M + 1], *FM = new double[M + 1];
	double *AN = new double[N + 1], *CN = new double[N + 1], *BN = new double[N + 1], *FN = new double[N + 1];
	double *lambd = new double[N + 1], *mu = new double[N + 1], *mas = new double[N + 1];
	const double R = param[0], r = param[1], dt = param[2] / T, a = param[3], dn = (2.0 * M_PI) / N, dm = (2.0 * M_PI) / M;
	const double ko = (dt*a*a) / (2.0*r*r), p = R / r;
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= M; j++)
		{
			u[0][i][j] = u0(dn*i, dm*j);
			data[0][j + i*(M + 1)] = u[0][i][j];
		}
		lambd[i] = sin(i*dn) / (p + cos(i*dn));
		mu[i] = 1.0 / ((p + cos(i*dn))*(p + cos(i*dn)));
	}
	for (int k = 0; k < T; k++)
	{
		for (int i = 0; i <= N; i++)
		{
			for (int j = 0; j <= M; j++)
			{
				fmas[k][i][j] = 0.5*dt*f(dn*i, dm*j, dt*k);
			}
		}
	}
	for (int k = 0; k<T - 1; k++)
	{
		for (int i = 0; i <= N; i++)
		{
			for (int j = 1; j<=M; j++)
			{
				AM[j] = -(mu[i] * ko) / (dm*dm);
				CM[j] = -1.0 + 2.0*AM[j];
				BM[j] = AM[j];
				FM[j] = -(u[k][i][j] + ko*((u[k][g(i + 1, N)][j] - 2 * u[k][i][j] + u[k][g(i - 1, N)][j]) / (dn*dn) - (lambd[i] / (2.0*dn)) * (u[k][g(i + 1, N)][j] - u[k][g(i - 1, N)][j])) + fmas[k][i][j]);
			}
			reduce(u_half[i], AM, CM, BM, FM, M);
		}
		for (int j = 0; j <= M; j++)
		{
			for (int i = 1; i <= N; i++)
			{
				AN[i] = (-ko / dn)*(1.0 / dn + 0.5*lambd[i]);
				CN[i] = -1.0 - (2.0*ko) / (dn*dn);
				BN[i] = (ko*lambd[i]) / (2.0*dn) - ko / (dn*dn);
				FN[i] = -(u_half[i][j] + ((mu[i] * ko) / (dm*dm))*(u_half[i][g(j + 1, M)] - 2.0*u_half[i][j] + u_half[i][g(j - 1, M)]) + fmas[k][i][j]);
			}
			reduce(mas, AN, CN, BN, FN, N);
			for (int i = 0; i <= N; i++)
				u[k + 1][i][j] = mas[i];
		}
		for (int i = 0; i <= N; i++)
		{
			for (int j = 0; j <= M; j++)
			{
				data[k + 1][j + i*(M + 1)] = u[k + 1][i][j];
			}
		}
	}
	delete[] mas;
	delete[] AN;
	delete[] CN;
	delete[] BN;
	delete[] FN;
	delete[] AM;
	delete[] CM;
	delete[] BM;
	delete[] FM;
	delete[] lambd;
	delete[] mu;
	for (int k = 0; k<T; k++)
	{
		for (int i = 0; i < N + 1; i++)
		{
			delete[] u[k][i];
			delete[] fmas[k][i];
		}
		delete[] u[k];
		delete[] fmas[k];
	}
	delete[] u;
	delete[] fmas;
	for (int j = 0; j < N + 1; j++)
		delete[] u_half[j];
	delete[] u_half;
}

