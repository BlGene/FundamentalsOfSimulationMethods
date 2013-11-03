#include <cstdio>
int main()
{
	float x = 0.01;
	double y = x;
	double z = 0.01;

	int i = x * 10000;
	int j = y * 10000;
	int k = z * 10000;

	printf("%d, %d, %d\n",i,j,k);

	return 0;
}
