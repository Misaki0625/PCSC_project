#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<math.h>
#include <cstdlib>
#include "function.h"
#include "readfile.h"
#include "Implement.h"
using namespace std;



Imple::Imple(const vector<float> v_m_read)
{
	this->v_read = v_m_read;
	int size = v_read.size();
	type = v_read[size-1]; // get the number of the last vector
	if (type == 1)
		cout << "The function type is polynomial" << endl;
	else if(type == 2)
		cout << "The function type is sinx & cosx" << endl;
	else
		cout << "The function type is exponential" << endl;
	v_read.pop_back();   // delete the last vector
	cout << "The read content is " << endl;
	printVector(v_read);
	v_input = function_input(v_read, type);
	int n = v_input.size();
	// Split v_input, one vector is placed on the left point, and the other vector is placed on the right point
	for (int i = 0; i < n / 2; i++)
	{
		v_input1.push_back(v_input[2 * i]);
		v_input2.push_back(v_input[2 * i + 1]);
	}
	cout << "the result for input is: " << endl;
	printVector(v_input1);
	printVector(v_input2);
	v_reinput = bisection_for_reinput(v_read, v_input1, v_input2, type);
	int m = v_reinput.size();
	for (int i = 0; i < m / 2; i++)
	{
		v_reinput1.push_back(v_reinput[2 * i]);
		v_reinput2.push_back(v_reinput[2 * i + 1]);
	}
	cout << "the result for reinput is: " << endl;
	printVector(v_reinput1);
	printVector(v_reinput2);
}

Equations::Equations(const Eigen::MatrixXf m)
{
	this->m_read = m;
	row = m_read.rows();
	Eigen::VectorXf v1(row);
	this->v_solution = v1;
	this->v_new = v1;
	for (int i = 0; i < row; i++)
	{
		v_solution[i] = 1;
	}
}

void Equations::func()
{

	int max_iter = 100;
	int iter = 0;
	float tol = 0.00001;
	for (int i = 0; i < max_iter; i++)
	{
		float diff = 0;
		Eigen::MatrixXf J(row, row);
		Eigen::VectorXf F(row);
		J = J_inverse(m_read, v_solution);
		F = f_non_linear_equations(m_read, v_solution);
		//printVector2(F);
		v_new = v_solution - J * F;
		for (int j = 0; j < row; j++)
		{
			diff += pow((v_new[j] - v_solution[j]), 2);
		}
		diff = sqrt(diff);
		v_solution = v_new;
		cout << "the diff is" << diff << endl;
		if (diff < tol)
		{
			cout << "The solution is " << endl;
			printVector2(v_new);
			cout << endl;
			cout << "The solve of the Function in this solution " << endl;
			printVector2(F);
			break;
		}
		float norm_f = 0;
		for (int j = 0; j < row; j++)
		{
			norm_f += pow(F[j], 2);
		}
		norm_f = sqrt(norm_f);
		//cout << norm_f << " * " << endl;
		//cout << "The F is " << endl;
		//printVector2(F);
		//cout << endl;
		if (norm_f < 0.01)
		{
			cout << "The solution is " << endl;
			printVector2(v_new);
			cout << endl;
			cout << "The solve of the Function in this solution " << endl;
			printVector2(F);
			break;
		}
		if (iter == (max_iter - 1))
		{
			cout << "not converge" << endl;
			printVector2(v_new);
			cout << endl;
			cout << "The solve of the Function in this solution " << endl;
			printVector2(F);
			break;
		}

		iter++;
	}
}

void Equations::get_output()
{
	output_EigenVector(v_new);
}

void Imple_bisection::func()
{
	this->v_out = bisection(v_read, v_reinput1, v_reinput2, type);
}

void Imple_bisection::get_output()
{
	printVector(v_out);
	output(v_out);
}

void Imple_chord::func()
{
	this->v_out = Classic_chord_method(v_read, v_reinput1, v_reinput2, type);
}

void Imple_chord::get_output()
{
	printVector(v_out);
	output(v_out);
}

void Imple_newton::func()
{
	this->v_out = Newton_method(v_read, v_reinput1, v_reinput2, type);
}

void Imple_newton::get_output()
{
	printVector(v_out);
	output(v_out);
}

void Imple_fixed_point::func()
{
	this->v_out = fixed_point_method(v_read, v_reinput1, v_reinput2, type);
}

void Imple_fixed_point::get_output()
{
	printVector(v_out);
	output(v_out);
}

void Imple_fixed_point_aitken::func()
{
	this->v_out = aitken_with_fixedpoint_method(v_read, v_reinput1, v_reinput2, type);
}

void Imple_fixed_point_aitken::get_output()
{
	printVector(v_out);
	output(v_out);
}