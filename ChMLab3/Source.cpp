#include <iostream>
#include <vector>
#include <string>

using namespace std;

void ShowProblem()
{
	cout << "д)" << endl;
	cout << "x^2/4 + y^2 - z = 0" << endl;
	cout << "x - z + 1 = 0" << endl;
}
vector<double> operator-(vector<double> a, vector<double> b)
{
	for (int i = 0; i < a.size(); ++i)
	{
		a[i] -= b[i];
	}
	return a;
}

vector<vector<double>> YakobiMatrix(double x, double y)
{
	vector<vector<double>> matrix;
	vector<double> temp;

	temp.push_back(x/2);
	temp.push_back(2 * y);
	matrix.push_back(temp);
	temp.clear();

	temp.push_back(1);
	temp.push_back(0);
	matrix.push_back(temp);

	return matrix;
}

double detMatrix(vector<vector<double>> a)
{
	return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

void transponation(vector<vector<double>>& a)
{
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = i + 1; j < a.size(); ++j)
		{
			double temp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = temp;
		}
	}
}

void Reverse(vector<vector<double>>& A)
{
	double detA = detMatrix(A);
	if ( detA == 0)
	{
		throw "Det A is 0";
	}
	double temp;
	transponation(A);
	for(int i = 0; i < A.size(); ++ i)
		for (int j = 0; j < A.size(); ++j)
		{
			A[i][j] /= detA;
		}
}

vector<double> matrixMultip(vector<vector<double>> a, vector<double> b)
{
	vector<double> c(b.size());
	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < b.size(); j++)
		{
			c[i] += a[i][j] * b[j];
		}
	}
	return c;
}

vector<vector<double>> matrixMultip(vector<vector<double>> a, vector<vector<double>> b)
{
	vector<vector<double>> c(a.size(), vector<double>(b[0].size()));
	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < b[0].size(); j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < a[0].size(); k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
	return c;
}

vector<double> Gauss(vector<vector<double>> a, vector<double> y, int n)
{
	double max;
	vector<double> x(n);
	int k, index;
	const double eps = 0.00001;
	k = 0;
	while (k < n)
	{
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		try
		{
			if (max < eps)
			{
				string error = "Нульовий рядок " + to_string(index) + " в матриці A\n";
				throw error;
			}
		}
		catch (string error)
		{
			cout << error;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;


		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue;
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue;
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}

vector<vector<double>> Invers(vector<vector<double>> a)
{
	vector<vector<double>> b;
	vector<double> e(a.size());
	for (int i = 0; i < a.size(); ++i)
	{
		e[i] = 1;
		b.push_back(Gauss(a, e, a.size()));
		e[i] = 0;
	}
	transponation(b);
	return b;
}



vector<double> F(double x0, double y0, double zm)
{
	double x1 = x0 * x0 / 4 + y0 * y0 - zm;
	double y1 = x0 - zm + 1;
	return { x1, y1 };
}

void showMatrix(vector<vector<double>> A)
{
	for (int i = 0; i < A.size(); ++i)
	{
		for (int j = 0; j < A.size(); ++j)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

double normVector(vector<double> a)
{
	double max = abs(a[0]);
	for (int i = 1; i < a.size(); ++i)
	{
		if (abs(a[i]) > max) max = abs(a[i]);
	}
	return max;
}
vector<double> Difference(vector<double> a, vector<double> b)
{
	vector<double> c;
	for(int i = 0; i < a.size(); ++i)
	{
		c.push_back(abs(a[i] - b[i]));
	}
	return c;
}

vector<double> modifNewton(double x0, double y0, double zm, double eps)
{
	vector<double> x0Vector = { x0, y0 }, x1Vector, difference;
	vector<vector<double>> A = YakobiMatrix(x0, y0);
	A = Invers(A);
	do
	{
		x1Vector = x0Vector - matrixMultip(A, F(x0Vector[0], x0Vector[1], zm));
		difference = Difference(x0Vector, x1Vector);
		x0Vector = x1Vector;
	} while (normVector(difference) > eps);
	return x1Vector;
}


int main()
{
	
	double x0 = 2, y0 = 2, zm = 1, eps = 1e-3, x;
	try
	{
		ShowProblem();
		for (zm = 1; zm < 6; ++zm)
		{
			vector<double> Ans = modifNewton(x0, y0, zm, eps);
			cout << "Ans for z = " << zm << " is  ";
			
			cout << "y = "<< Ans[1] << " " ;
			cout << "and x = " << Ans[0];
			cout << endl << endl;
			x0++;
			y0++;
		}
		
	}
	catch (...)
	{
		cout << "Det A is 0";
	}
	
}