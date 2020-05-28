#include <iostream>
#include <math.h>
#include <fstream>
#include<string>
using namespace std;

void coefC(double* cArr, double** mass, int n, double h)
{
	int i, j;
	double** A = new double* [n - 2];
	for (i = 0; i < n - 2; i++)
		A[i] = new double[n - 2];

	for (i = 1; i < n - 3; i++)
	{
		j = i - 1;
		A[i][j] = h;
		A[i][j + 1] = 4 * h;
		A[i][j + 2] = h;
	}
	A[0][0] = 4 * h;
	A[0][1] = h;
	A[n - 3][n - 4] = h;
	A[n - 3][n - 3] = 4 * h;

	double** F = new double* [n - 2];
	for (i = 0; i < n - 2; i++)
		F[i] = new double[1];
	

	for (i = 1; i < n - 1; i++)
		F[i - 1][0] = 3 * ((mass[i + 1][1] - mass[i][1]) / h - (mass[i][1] - mass[i - 1][1]) / h);

	for (i = 1; i < n - 2; i++)
	{
		F[i][0] -= (A[i][i - 1] * F[i - 1][0] / A[i - 1][i - 1]);
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1]);
	}

	for (i = 1; i < n - 2; i++)
	{
		F[i][0] -= (A[i][i - 1] * F[i - 1][0] / A[i - 1][i - 1]);
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1]);
	}

	cArr[n - 2] = F[n - 3][0] / A[n - 3][n - 3];
	for (i = n - 4; i >= 0; i--)
		cArr[i + 1] = (F[i][0] - A[i][i + 1] * cArr[i + 2]) / A[i][i];
	
	cArr[0] = 0;
	cArr[n - 1] = 0;
	for (int i = 0; i < n - 2; i++)
	{
		delete[] A[i];
		delete[] F[i];
	}
	delete[] A;
	delete[] F;
}

//spline
double spline(double x, double* a, double* b, double* p, double* q, double** func, int n)
{
	int i = 0;
	for (int j = 1; j < n - 1; j++)
	{
		if (x > func[j][0])
			i++;
		else
			break;
	}
	return a[i] + b[i] * (x - func[i][0]) + p[i] * pow((x - func[i][0]), 2) + q[i] * pow((x - func[i][0]), 3);
}

int main()
{
	int i, n = 9;
	double** f;
	f = new double* [n];
	double* a = new double[n - 1];
	double* b = new double[n - 1];
	double* c = new double[n];
	double* d = new double[n - 1];

	for (int i = 0; i <= n - 1; i++)
	{
		f[i] = new double[n];
		for (int j = 0; j <= 1; j++)
		{
			if ((i == 0) && (j == 0))
				f[i][j] = 0;
			else if ((i == 0) && (j == 1))
				f[i][j] = 0;
			else if ((i == 1) && (j == 0))
				f[i][j] = 0.125;
			else if ((i == 1) && (j == 1))
				f[i][j] = 0.12467;
			else if ((i == 2) && (j == 0))
				f[i][j] = 0.25;
			else if ((i == 2) && (j == 1))
				f[i][j] = 0.247234;
			else if ((i == 3) && (j == 0))
				f[i][j] = 0.375;
			else if ((i == 3) && (j == 1))
				f[i][j] = 0.364902;
			else if ((i == 4) && (j == 0))
				f[i][j] = 0.5;
			else if ((i == 4) && (j == 1))
				f[i][j] = 0.473112;
			else if ((i == 5) && (j == 0))
				f[i][j] = 0.625;
			else if ((i == 5) && (j == 1))
				f[i][j] = 0.563209;
			else if ((i == 6) && (j == 0))
				f[i][j] = 0.75;
			else if ((i == 6) && (j == 1))
				f[i][j] = 0.616193;
			else if ((i == 7) && (j == 0))
				f[i][j] = 0.875;
			else if ((i == 7) && (j == 1))
				f[i][j] = 0.579699;
			else if ((i == 8) && (j == 0))
				f[i][j] = 1;
			else if ((i == 8) && (j == 1))
				f[i][j] = 0;
		}
	}

	
	double h = 0.125;
	coefC(c, f, n, h);
	for (i = 0; i < n - 1; i++)
	{
		a[i] = f[i][1];
		b[i] = (f[i + 1][1] - f[i][1]) / h - c[i] * h - (c[i + 1] - c[i]) * h / 3.;
		d[i] = (c[i + 1] - c[i]) / (3 * h);
	}

	//method trapeciy
	double up = 0, down = 1;
	int ndrugaai = 999999;
	double width = (down - up) / ndrugaai;

	double trapez_integr = 0;
	for (int step = 0; step < ndrugaai; step++) {
		const double x1 = up + step * width;
		const double x2 = up + (step + 1) * width;

		trapez_integr += 0.5 * (x2 - x1) * (spline(x1, a, b, c, d, f, n) + spline(x2, a, b, c, d, f, n));
	}

	//simpson method
	ndrugaai = 999999;
	width = (down - up) / ndrugaai;
	double simps_integr = 0;
	for (int step = 0; step < ndrugaai; step++) {
		const double x1 = up + step * width;
		const double x2 = up + (step + 1) * width;
		simps_integr += (x2 - x1) / 6.0 * (spline(x1, a, b, c, d, f, n) + 4.0 * spline(0.5 * (x1 + x2), a, b, c, d, f, n) + spline(x2, a, b, c, d, f, n));
	}

	//runge rule for trapecii method
	double trapec;
	double tochnost_trapecii;
	double polovina_trapeciya;
	double ndrugaaitrap;
	trapec = trapez_integr;
	ndrugaai = 999999;
	ndrugaaitrap = ndrugaai / 2;
	width = (down - up) / ndrugaaitrap;
	trapez_integr = 0;
	for (int step = 0; step < ndrugaaitrap; step++)
	{
		const double x1 = up + step * width;
		const double x2 = up + (step + 1) * width;
		trapez_integr += 0.5 * (x2 - x1) * (spline(x1, a, b, c, d, f, n) + spline(x2, a, b, c, d, f, n));
	}
	polovina_trapeciya = trapez_integr;
	tochnost_trapecii = fabs((trapec- polovina_trapeciya) / trapec);
	cout << "metod trapecii: " << trapez_integr << "\t" << "pravilo rynge dlai metoda trapecii: " << tochnost_trapecii << "\t" << endl;

	//ans1.dat
	ofstream trap("ans1.dat");
	trap << trapec << "\t" << tochnost_trapecii << endl;
	trap.close();

	//runge rule for simpson method
	double simpson;
	double tochnost_simpsona;
	double polovina_simpsona;
	double ndrugaaisimp;
	simpson = simps_integr;
	ndrugaai = 999999;
	ndrugaaisimp = ndrugaai / 2;
	width = (down - up) / ndrugaaisimp;

	simps_integr = 0;
	for (int step = 0; step < ndrugaaisimp; step++) {
		const double x1 = up + step * width;
		const double x2 = up + (step + 1) * width;

		simps_integr += (x2 - x1) / 6.0 * (spline(x1, a, b, c, d, f, n) + 4.0 * spline(0.5 * (x1 + x2), a, b, c, d, f, n) + spline(x2, a, b, c, d, f, n));
	}
	polovina_simpsona = simps_integr;

	tochnost_simpsona = fabs((simpson - polovina_simpsona) / simpson);

	cout << "metod simpsona: " << simpson << "\t" << "pravilo rynge dlai metoda simpsona: " << tochnost_simpsona << "\t" << endl;

	//ans2.dat
	ofstream simp("ans2.dat");
	simp << simpson << "\t" << tochnost_simpsona << endl;
	simp.close();
	return 0;
}