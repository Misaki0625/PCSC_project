#include <iostream>
#include <vector>
#include "function.h"
#include "implement.h"
#include "math.h"
#include <stdlib.h>  

using namespace std;

float chord_test(string path)
{
	float diff = 0;
	vector<float> v_read = Read_File_txt(path);
	Imple_chord r1(v_read);
	r1.func();



	return diff;
}