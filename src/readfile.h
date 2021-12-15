#pragma once
#ifndef __READFILE__H__
#define __READFILE__H__
#include <iostream>
#include <Eigen/Dense>
using namespace std;


/**
* This is my ReadFile Class
*/
class ReadFile
{
/// <summary>
/// This class is used to the be Pure virtual class for reading files
/// </summary>
public:

	virtual vector<float> read() = 0;

	string file_path;
};


/**
* This is my ReadFile_txt Class
*/
class ReadFile_txt : private ReadFile
{
/// <summary>
///  The ReadFile_txt class is used to read polynomial covariants from txt files. If we read the number a, b, c, d, it means we get a polynomial function of ax ^ 3 + bx ^ 2 + cx + d
/// </summary>
public:

	ReadFile_txt(const string path);

	virtual vector<float> read();

};


/**
* This is my ReadFile_csv Class
*/
class ReadFile_csv : private ReadFile
{
/// <summary>
/// The ReadFile_csv class is used to read polynomial covariants from txt files. If we read the number a, b, c, d,e it means we get a polynomial function of ax ^ 4 + bx ^ 3 + cx^2 + dx + e
/// </summary>
public:
	ReadFile_csv(const string path);

	virtual vector<float> read();

};

/**
* This is my ReadFile_matrix_txt Class
*/
class ReadFile_matrix_txt : private ReadFile
{
	/// <summary>
	///  The ReadFile_matrix_txt class is used to read non linear equations covariants from txt files. 
	/// </summary>
public:

	virtual vector<float> read();

	ReadFile_matrix_txt(const string path);

	virtual Eigen::MatrixXf read_matrix();

};
#endif