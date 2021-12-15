#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "readfile.h"
#include <Eigen/Dense>

using namespace std;


ReadFile_txt::ReadFile_txt(const string path)
{
	this->file_path = path;
}


vector<float> ReadFile_txt::read()
{
	vector<float> s;
	vector<float> k;
	k.push_back(0);
	ifstream ifs;
	ifs.open(file_path, ios::in);
	if (!ifs.is_open())
	{
		cout << "failed to open the file: " << file_path << endl;
		return k; //If we fail to open the file, it will return a vector contains zero
	}
	for (float a; ifs >> a;)
	{
		s.push_back(a); // We put the number we read to a vector
	}
	ifs.close();

	return s;
	}


ReadFile_csv::ReadFile_csv(const string path)
{
	this->file_path = path;
}


vector<float> ReadFile_csv::read()
{
	vector<float> s;
	vector<float> k;
	k.push_back(0);
	vector<string> row;
	string line, word;
	string fname = "test01.csv";
	vector<vector<string>> content;
	fstream file(fname, ios::in);
	if (!file.is_open())
	{
		cout << "failed to open the file: " << file_path << endl;
		return k; //If we fail to open the file, it will return a vector contains zero
	}
	if (file.is_open())
	{
		while (getline(file, line))
		{
			row.clear();

			stringstream str(line);

			while (getline(str, word, ','))
				row.push_back(word);
			content.push_back(row);
		}
	}
	else
		cout << "Could not open the file\n";

	for (int i = 0; i < content.size(); i++)
	{
		for (int j = 0; j < content[i].size(); j++)
		{
			s.push_back(std::stod(content[i][j]));
		}
	}
	file.close();
	return s;
}

Eigen::MatrixXf ReadFile_matrix_txt::read_matrix()
{
	Eigen::MatrixXf k(1,1);
	k(0,0) = 0;
	vector<float> v;

	ifstream ifs;
	ifs.open(file_path, ios::in);
	if (!ifs.is_open())
	{
		cout << "failed to open the file: " << file_path << endl;
		return k; //If we fail to open the file, it will return a vector contains zero
	}
	for (float a; ifs >> a;)
	{
		v.push_back(a); // We put the number we read to a vector
		//cout << a<<" ";
	}
	int len = v.size();
	//cout << len << endl;
	int n = (-1 + sqrt(1 + 12 * len)) / 6; // rows

	Eigen::MatrixXf m(n, 3 * n + 1);
	//cout << m.size() << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < (3 * n + 1); j++)
		{
			m(i, j) = v[i * (3 * n + 1) + j];
		}
	}
	ifs.close();
	return m;
}

ReadFile_matrix_txt::ReadFile_matrix_txt(const string path)
{
	this->file_path = path;
}

vector<float> ReadFile_matrix_txt::read()
{
	vector<float> v;
	v.push_back(0);
	cout << "please use the read_matrix function" << endl;
	return v;
}