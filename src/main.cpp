#include <iostream>
#include <vector>
#include "function.h"
#include "implement.h"
#include "math.h"
#include <stdlib.h>  

using namespace std;


int main() // solve the zero point for polynomial function with different implementations
{
	// path loading
	string path_poly = "test.txt";
	string path_poly_csv = "test01.csv";
	string path_traingular = "test_tri.txt";
	string path_exponential = "test_exp.txt";
	string path_matrix = "matrix.txt";
	
	//reading part
	//vector<float> v_read = Read_File_txt(path_traingular); // choose between reading a txt or csv file, and choose the path you want to read(except for path_matrix)
	//vector<float> v_read = Read_File_csv(path_poly);
	Eigen::MatrixXf m_read = Read_matrix_txt(path_matrix); // used to read the path_matrix

	// Implementation part
	//Imple_chord r1(v_read);
	//Imple_bisection r1(v_read);
	//Imple_newton r1(v_read);
	//Imple_fixed_point r1(v_read);
	//Imple_fixed_point_aitken r1(v_read);
	Equations r1(m_read); // only used for the path_matrix reading
	
	// execute part
	r1.func();
	r1.get_output();

	return 0;
}
