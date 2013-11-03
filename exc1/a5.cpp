#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <gmpxx.h>

using namespace std;

//comparison to sort values by magnitude
bool compare(double a, double b)
{
	return (fabs(a)<fabs(b));
}

int main()
{
	//stream to read in the numbers
	ifstream is ("numbers.dat", ifstream::binary);
	//determine size (in bytes)
	is.seekg (0, is.end);
	int length = is.tellg();
	is.seekg (0, is.beg);
	//allocate memory to read in the file
	char * buffer = new char [length];

	cout << "Reading " << length << " characters... ";
	// read data as a block:
	is.read (buffer,length);
	if (is)
		cout << "all characters read successfully." << endl;
	else
		cout << "error: only " << is.gcount() << " could be read" << endl;
	//closing the stream
	is.close();
	
	//printing the first integer = #values in file
	int* file_len_p= reinterpret_cast<int*>(buffer);
	cout << "There are " << *file_len_p << " doubles in this file." << endl;

	//allocating memory to cast character array read in to double array
	double * dbuff = new double[*file_len_p];
	//casting character array to double array via memcopy
	memcpy(dbuff,buffer+sizeof(int),*file_len_p*sizeof(double));
	//creating a vector from this array via assign for convenience
	vector<double> vec;
	vec.assign(dbuff,dbuff+*file_len_p);

	cout << "Summing up all the values according to excercise a) - e) with results:" << endl; 
	//a) summing up the values from beginning to end
	double sum_a = 0.;
	for(auto x : vec)
		sum_a += x;
	cout << "Part a): " << sum_a << endl;

	//b) summing up the values from end to beginning
	double sum_b = 0.;
	for(auto j = vec.end(); j != vec.begin(); j--)
		sum_b += *j;
	cout << "Part b): " << sum_b << endl;

	//c) sorting the values by magnitude and then summing up from small to big
	//copying the vector and then sorting the entries via std::sort
	auto vec_sort = vec;
	sort(vec_sort.begin(),vec_sort.end(),compare);
	double sum_c = 0.;
	for(auto x : vec_sort)
		sum_c += x;	
	cout << "Part c): " << sum_c << endl;

	//d) summing up the sorted values with a long double 
	long double sum_d = 0.;
	for(auto x : vec_sort)
		sum_d += x;	
	cout << "Part d): " << sum_d << endl;

	//e) summing up the sorted values with the gmp library
	mpf_t sum_e;
	mpf_t summand;
	mpf_init2(sum_e,10);
	mpf_init2(summand,10);
	mpf_set_d(sum_e,0.);
	for(auto x : vec_sort)
	{
		mpf_set_d(summand,x);
		mpf_add(sum_e,sum_e,summand);
	}
	gmp_printf("Part e): %. *Ff", 5, sum_e);
	cout << endl;

	return 0;
}
