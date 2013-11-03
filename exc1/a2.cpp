#include <iostream>

int main()
{
	double a = 1.0e17;
	double b = -1.0e17;
	double c = 1.0;

	double x = (a + b) + c;
	std::cout << "x = " << x << std::endl;

	double y  = a + (b + c);
	std::cout << "y = " << y << std::endl;

	return 0;
}
