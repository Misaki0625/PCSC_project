#pragma once
#ifndef __FUNCTION__H__
#define __FUNCTION__H__
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<math.h>
#include <cstdlib>
#include <Eigen/Dense>

using namespace std;


float f(vector<float> v, int n, float x, float type);

float derivative_f(vector<float> v, int n, float x, float type);

Eigen::VectorXf f_non_linear_equations(Eigen::MatrixXf m, Eigen::VectorXf vec);

Eigen::MatrixXf J_inverse(Eigen::MatrixXf m, Eigen::VectorXf vec);

void myswap(float* a, float* b);

void printVector(vector<float>& v);

void printVector2(Eigen::VectorXf vec);

void printMatrix(Eigen::MatrixXf m);

vector<float> Read_File_txt(string file_path);

vector<float> Read_File_csv(string file_path);

Eigen::MatrixXf Read_matrix_txt(string file_path);

void output(vector<float>& v);

void output_EigenVector(Eigen::VectorXf v);

vector<float> function_input(vector<float>& v, float type);

vector<float> bisection_for_reinput(vector<float>& v, vector<float>& v1, vector<float>& v2, float type);

vector<float> bisection(vector<float>& v, vector<float>& v1, vector<float>& v2, float type);

vector<float> Classic_chord_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type);

vector<float> Newton_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type);

vector<float> fixed_point_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type);

vector<float> aitken_with_fixedpoint_method(vector<float>& v, vector<float>& v1, vector<float>& v2, float type);  // Steffensen


#endif