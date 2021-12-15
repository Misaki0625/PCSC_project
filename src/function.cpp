#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<math.h>
#include <cstdlib>
#include "function.h"
#include "readfile.h"
using namespace std;


/*!
\param v an vector argument.
\param n a int argument which is the size of v.
\param x a float argument 
\return a, the result of f(x), for which the v is the covariants 
f function description
*/
float f(vector<float> v, int n, float x, float type) { 
	
	// To get the number of f(x0) with the x0 and the function coefficient provided
	if (type == 1)
	{
		float a = 0;
		for (int i = 0; i < n; i++)
		{
			a += v[i] * pow(x, n - i - 1);
		}
		return a;
	}
	else if (type == 2)
	{
		float a = 0;
		a = v[0] * sin(v[1] * x + v[2]) + v[3] * cos(v[4] * x + v[5]) + v[6];
		return a;
	}
	else
	{
		float a = 0;
		for (int i = 0; i < (n - 1) / 2; i++)
		{
			a += v[2*i] * exp(v[2*i+1]*x);
		}
		a += v[n-1];
		return a;
	}
	
}


/*!
\param v an vector argument.
\param n a int argument which is the size of v.
\param x a float argument
\return a, the result of f'(x), for which the v is the covariants
derivative_f function description
*/
float derivative_f(vector<float> v, int n, float x, float type)
{
	if (type == 1)
	{
		// To get the number of the derivative of f(x0) with the x0 and the function coefficient provided
		float a = 0;
		for (int i = 0; i < n - 1; i++)
		{
			a += v[i] * pow(x, n - i - 2) * (n - i - 1);
		}
		return a;
	}
	else if (type == 2)
	{
		float a = 0;
		a = -v[3] * v[4] * sin(v[4] * x + v[5]) + v[0] * v[1] * cos(v[1] * x + v[2]) + v[6];
		return a;
	}
	else
	{
		float a = 0;
		for (int i = 0; i < (n - 1) / 2; i++)
		{
			a += v[2*i] * v[2*i+1] * exp(v[2*i + 1] * x);
		}
		return a;
	}


}

Eigen::VectorXf f_non_linear_equations(Eigen::MatrixXf m, Eigen::VectorXf vec)
{
	// v is the solution of x1,x2,...,xn
	int row = m.rows();
	Eigen::VectorXf v_solve(row);
	for (int i = 0; i < row; i++)
	{
		float a = 0;
		for (int j = 0; j < row; j++)
		{
			a += m(i, 3 * j) * pow(vec[j], 3) + m(i, 3 * j + 1) * pow(vec[j], 2) + m(i, 3 * j + 2) * vec[j]; // 3a*x^3 + 2b*x^2 + c*x
		}
		a += m(i, 3 * row);
		v_solve[i] = a;
	}
	return v_solve;
}

Eigen::MatrixXf J_inverse(Eigen::MatrixXf m, Eigen::VectorXf vec)
{
	int col = m.cols();
	int row = m.rows();
	Eigen::MatrixXf J(row, row); // shape : n * n
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < row; j++)
		{
			J(i, j) = 3 * m(i, 3 * j) * pow(vec[j], 2) + 2 * m(i, 3 * j + 1) * vec[j] + m(i, 3 * j + 2); // 3a*x^2 + 2b*x + c
		}

	}
	J = J.inverse();

	return J;

}


void myswap(float* a, float* b)
{
	float t = *a;
	*a = *b;
	*b = t;
}


/*!
\param v an vector argument.
printVector function description
*/
void printVector(vector<float>& v) { 
	
	// Print the component in the vector
	for (vector<float>::iterator it = v.begin(); it != v.end(); it++) {
		cout << *it << " ";
	}
	cout << endl;
}

void printVector2(Eigen::VectorXf vec)
{
	int row = vec.rows();
	for (int k = 0; k < row; k++)
	{
		cout << vec[k] << " ";
	}
}

void printMatrix(Eigen::MatrixXf m)
{
	int row = m.rows();
	int col = m.cols();
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			cout << m(i, j)<<" ";
		}
		cout << endl << endl;;
	}
}

/*!
\param file_path, the URL we need to read the vector
\return v0, a set of number for coefficients for polynomial function
Read_File function description
*/
vector<float> Read_File_txt(string file_path)
{
	vector<float> v_read;
	ReadFile_txt r1 = ReadFile_txt(file_path);
	v_read = r1.read();
	return v_read;
}

/*!
\param file_path, the URL we need to read the vector
\return v0, a set of number for coefficients for polynomial function
Read_File function description
*/
vector<float> Read_File_csv(string file_path)
{
	vector<float> v_read;
	ReadFile_csv r1 = ReadFile_csv(file_path);
	v_read = r1.read();
	return v_read;
}


Eigen::MatrixXf Read_matrix_txt(string file_path)
{
	Eigen::MatrixXf m_read;
	ReadFile_matrix_txt r1 = ReadFile_matrix_txt(file_path);
	m_read = r1.read_matrix();
	cout << "The read matrix is" << endl;
	printMatrix(m_read);
	return m_read;
}

/*!
\param v an vector argument.
write the component in the vector to a txt file
output_csv function description
*/
void output(vector<float>& v)
{
	// Write the zero points we get to the txt file
	ofstream ofs;
	ofs.open("output.txt", ios::out);
	//ofs <<"The solution(zero point) of the function is :" << endl;
	for (vector<float>::iterator it = v.begin(); it != v.end(); it++) {
		ofs << *it << " ";
	}
	cout << endl;
	ofs.close();
}


void output_EigenVector(Eigen::VectorXf v)
{
	// Write the zero points we get to the txt file
	int n = v.rows();
	ofstream ofs;
	ofs.open("output.txt", ios::out);
	//ofs <<"The solution(zero point) of the function is :" << endl;
	for (int i = 0; i < n;i++) {
		ofs << v[i] << " ";
	}
	cout << endl;
	ofs.close();
}


/*!
\param v an vector argument.
\return v1, the end point of the interval where the zero point lies
function_input function description
*/
vector<float> function_input(vector<float>& v, float type) 
{
	vector<float> v1; // Storage point, even index is left point, odd index is right point
	int n = v.size();
	if (type == 2)
	{
		float a = 2 * 3.1416 / v[1];
		float b = 2 * 3.1416 / v[4];
		float c = a * b;
		for (float i = 0; i < c + 1; i = i + c/100) //we define the zero point search area in [-1000000,100000]
		{
			if (f(v, n, i, type) * f(v, n, i + 0.1, type) < 0) {
				v1.push_back(i);
				v1.push_back(i + 0.1);
			}
		}
		return v1;
	}
	else
	{
		for (float i = -100000; i < 100000; i = i + 0.1) //we define the zero point search area in [-1000000,100000]
		{
			if (f(v, n, i, type) * f(v, n, i + 0.1, type) < 0) {
				v1.push_back(i);
				v1.push_back(i + 0.1);
			}
		}
		return v1;
	}

}


/*!
\param v an vector argument, the covariants of the polynomial function.
\param v1 an vector argument, a set of left point for the intervals where the zero point lies.
\param v2 an vector argument, a set of right point for the intervals where the zero point lies
\return v0, a set of left and right point for the minimized interval
bisection_for_reinput function description
*/
vector<float> bisection_for_reinput(vector<float>& v, vector<float>& v1, vector<float>& v2, float type)
{

	// To minimize the interval for zero point from 0.1 to 0.001
	vector<float> v0; //Store the left and right point for the minimized interval we get. Even index is left point, odd index is right point.

	int n = v.size();

	int length = v1.size();

	float mid = 0;

	if (length == 0)
	{
		cout << "There may be no solution between [-100,000, 100,000)" << endl << endl;
	}
	for (int i = 0; i < length; i++) 
	{
		
		float x1 = v1[i];
		float x2 = v2[i];
		// cout << x1 << " " << x2 << endl;
		while (abs(float(x2 - x1)) > 0.001) // minimize the interval to 0.001
		{
			mid = (x1 + x2) / 2;
			if (f(v, n, mid, type) * f(v, n, x2, type) < 0)
				x1 = mid;
			else
				x2 = mid;
			// cout << x1 << " " << x2 << endl;
		}

		v0.push_back(x1);
		v0.push_back(x2);

	}

	return v0;
}


/*!
\param v an vector argument, the covariants of the polynomial function.
\param v1 an vector argument, a set of left point for the intervals where the zero point lies.
\param v2 an vector argument, a set of right point for the intervals where the zero point lies
\return v0, a set of zero points
bisection function description
*/
vector<float> bisection(vector<float>& v, vector<float>& v1, vector<float>& v2, float type) // bisection求f(x)=0的零点
{
	vector<float> v0;

	int n = v.size();
    
	int length =v1.size();

	float mid = 0;

	for (int i = 0; i < length; i++) // use the bisection to get the zeros point
	{
		
		float x1 = v1[i];
		float x2 = v2[i];
		//cout << x1 << " " << x2 << endl;
		while (abs(float(x2 - x1)) > 0.000001) // use the bisection to get the zeros point
		{ 
			
			mid = (x1 + x2) / 2;
			if (f(v, n, mid, type) * f(v, n, x2, type) < 0)
				x1 = mid;
			else
				x2 = mid;
		}
		
		v0.push_back(mid);
		
	}
	cout << "The point solved by the bisection is : " << endl;
	return v0;
}


/*!
\param v an vector argument, the covariants of the polynomial function.
\param v1 an vector argument, a set of left point for the intervals where the zero point lies.
\param v2 an vector argument, a set of right point for the intervals where the zero point lies
\return v0, a set of zero points
Classic_chord_method function description
*/
vector<float> Classic_chord_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type) // 弦法
{
	vector<float> v0;
	int n = v.size();
	int length = v1.size();
	float x = 0;
	for (int i = 0; i < length; i++)
	{
		float x1 = v1[i];
		float x2 = v2[i];
		float a = v1[i];
		float b = v2[i];
		while (abs(float(x2 - x1)) > 0.0001)
		{
			x = x2 - f(v, n, x2, type) * (x2 - x1) / (f(v, n, x2, type) - f(v, n, x1, type)) ; // classic chord function for update
			if (f(v, n, x2, type) * f(v, n, x, type) < 0)
			{
				x1 = x;
			}
			else
			{
				x2 = x;
			}
		}
		if(x >= a && x <= b)
		{
			v0.push_back(x);
		}
		else
		{
			cout << "The function divergent in  interval [" << a << "," << b << "]" << endl;
		}
	}
	cout << "The point solved by the classic chord method is : " << endl;
	return v0;
}


/*!
\param v an vector argument, the covariants of the polynomial function.
\param v1 an vector argument, a set of left point for the intervals where the zero point lies.
\param v2 an vector argument, a set of right point for the intervals where the zero point lies
\return v0, a set of zero points
Newton_method function description
*/
vector<float> Newton_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type)
{
	vector<float> v0;
	int n = v.size();
	int length = v1.size();
	for (int i = 0; i < length; i++)
	{
		float x = (v1[i] + v2[i])/2;
		float a = v1[i];
		float b = v2[i];
		int m = 0;
		while (abs(f(v, n, x, type)) < 0.0001)
		{
			x = x - f(v, n, x, type)/ derivative_f(v, n, x, type); // Newton–Raphson method for update			
			m++;
			if (x < a || x > b) // if x is in the interval
			{
				cout << "The function divergent in  interval ["<<a<<","<<b<<"]" << endl;
				x = b + 0.01;
				break;
			}
			if (m == 100)
			{
				cout << "The function divergent in  interval [" << a << "," << b << "]" << endl;
				
				break;
			}
		}
		v0.push_back(x);
	}
	cout << "The point solved by the Newton–Raphson method is : " << endl;
	return v0;
}


/*!
\param v an vector argument, the covariants of the polynomial function.
\param v1 an vector argument, a set of left point for the intervals where the zero point lies.
\param v2 an vector argument, a set of right point for the intervals where the zero point lies
\return v0, a set of zero points
fixed_point_method function description
*/
vector<float> fixed_point_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type)
{
	vector<float> v0;
	int n = v.size();
	int length = v1.size();
	float x;
	for (int i = 0; i < length; i++)
	{
		float x = (v1[i]+v2[i])/2;
		float a = v1[i];
		float b = v2[i];
		float x1 = x + 1; // it's only enable that we can step into the first circulation
		while (abs(float(x - x1)) > 0.0001)
		{
			x1 = x;
			x = f(v, n, x1, type) + x1; //g(x) = f(x) - x  ;   x' = g(x) -> x''= g(x')
			if (abs(x - x1) > 1)// use to find whether if it diverges
			{
				cout << "The function divergent in  interval [" << a << "," << b << "]" << endl;
				break;
			}

		}
		if (x >= a && x <= b)
		{
		v0.push_back(x);
		}
	}
	cout << "The point solved by the fixed point method is : " << endl;
	return v0;
}


/*!
\param v an vector argument, the covariants of the polynomial function.
\param v1 an vector argument, a set of left point for the intervals where the zero point lies.
\param v2 an vector argument, a set of right point for the intervals where the zero point lies
\return v0, a set of zero points
aitken_with_fixedpoint_method function description
*/
vector<float> aitken_with_fixedpoint_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type) 
{	// It is also caleed Steffensen method
	// => g(x) = x ; f(x) = 0 ; g(x) = f(x) + x
	// y_k = g(x_k)
	// z_k = g(y)
	// x_k+1 = x_k - (y_k - x_k)^2/(z_k - 2*y_k +x_k)
	
	vector<float> v0;
	int n = v.size();
	int length = v1.size();
	float x;
	for (int i = 0; i < length; i++)
	{
		float a = v1[i];
		float b = v2[i];
		x = (v1[i] + v2[i]) / 2;
		float y = 0;
		float z = 0;
		float x1 = x + 1; 
		while (abs(float(x - x1)) > 0.000001)
		{
			x1 = x;
			y = f(v, n, x1, type) + x1;
			z = f(v, n, y, type) + y;
			x = x1 - (y - x1) * (y - x1) / (z - 2 * y + x);
			if (abs(x - x1) > 1)// use to find whether if it diverges
			{
				cout << "The function divergent in  interval [" << a << "," << b << "]" << endl;
				break;
			}
		}

		if (x >= a && x <= b)
		{
			v0.push_back(x);
		}
	}
	cout << "The point solved by the fixed point method with aitken acceleration is : " << endl;
	return v0;
}