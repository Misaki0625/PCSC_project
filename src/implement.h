#pragma once
#ifndef __IMPLEMENT__H__
#define __IMPLEMENT__H__
#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace std;


/**
* This is my Imple Class
*/
class Imple
{
	/// <summary>
	/// This class is used to do some preprocess on the function
	/// </summary>

public:
	Imple(const vector<float> v_m_read);

	vector<float> v_read, v_input, v_input1, v_input2, v_reinput, v_reinput1, v_reinput2;

	float type;

};


class Equations
{
public:
	Equations(const Eigen::MatrixXf m);

	void func();

	void get_output();

private:
	int row;
	Eigen::MatrixXf m_read;
	Eigen::VectorXf v_solution;
	Eigen::VectorXf v_new;

};

/**
* This is my Imple_bisection Class
*/
class Imple_bisection : public Imple
{
	/// <summary>
	/// This class is used to use the bisection method to figure out the zero point
	/// </summary>
public:

	Imple_bisection(vector<float> v_m_read) : Imple(v_m_read)
	{
		this->v_read = v_m_read;
	}

	void func(); // run the implementation

	void get_output();
private:
	vector<float> v_out;
};

/**
* This is my Imple_chord Class
*/
class Imple_chord : public Imple
{
	/// <summary>
	/// This class is used to use the classic chord method to figure out the zero point
	/// </summary>
public:

	Imple_chord(vector<float> v_m_read) : Imple(v_m_read)
	{
		this->v_read = v_m_read;
	}

	void func();

	void get_output();
private:
	vector<float> v_out;
};


/**
* This is my Imple_newton Class
*/
class Imple_newton : public Imple
{
	/// <summary>
	/// This class is used to use the newton method to figure out the zero point
	/// </summary>
public:

	Imple_newton(vector<float> v_m_read) : Imple(v_m_read)
	{
		this->v_read = v_m_read;
	}

	void func();

	void get_output();
private:
	vector<float> v_out;
};


/**
* This is my Imple_fixed_point Class
*/
class Imple_fixed_point : public Imple
{
	/// <summary>
	/// This class is used to use the fixed point method to figure out the zero point
	/// </summary>
public:

	Imple_fixed_point(vector<float> v_m_read) : Imple(v_m_read)
	{
		this->v_read = v_m_read;
	}

	void func();

	void get_output();
private:
	vector<float> v_out;
};


/**
* This is my Imple_fixed_point_aitken Class
*/
class Imple_fixed_point_aitken : public Imple
{
	/// <summary>
	/// This class is used to use the fixed point method with aitken acceleration to figure out the zero point
	/// </summary>
public:

	Imple_fixed_point_aitken(vector<float> v_m_read) : Imple(v_m_read)
	{
		this->v_read = v_m_read;
	}

	void func();

	void get_output();
private:
	vector<float> v_out;
};

#endif

