#include <iostream>

//templated function to determine the machine precision
template<class T>
T get_precision(T starting_value, T divisor)
{
	//initializing eps with the starting value (typically 1)
	T eps = starting_value;
	//divide eps by divisor as long starting_value + eps yields a result not equal to starting value
	while(starting_value + eps != starting_value)
	{
		eps /= divisor;
	}	
	//multiply eps with the divisor because the value after the loop is below machine precision
	eps *= divisor;
	return eps;
}

//for the divisor a value near to 1 is chosen to obtain a value as near to the "real eps" as possible
//the machine precision is determined for float, double and long double
int main()
{
	std::cout << "Determining epsilon of float, which has " << sizeof(float) << " bytes" << std::endl;
	std::cout << "Result: eps_float = " << get_precision<float>(1.,1.000001) << std::endl;

	std::cout << "Determining epsilon of double, which has " << sizeof(double) << " bytes" << std::endl;
	std::cout << "Result: eps_double = " << get_precision<double>(1.,1.000001) << std::endl;
	
	std::cout << "Determining epsilon of long double, which has " << sizeof(long double) << " bytes" << std::endl;
	std::cout << "Result: eps_long = " << get_precision<long double>(1.,1.000001) << std::endl;
	
	return 0;
}
