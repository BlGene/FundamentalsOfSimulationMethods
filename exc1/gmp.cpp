#include <gmp.h>


int main(int rgc, char *argv[])
{
	mpf_t a, b, c;

	mpf_init2(a, 30);
	mpf_init2(b, 30);
	mpf_init2(c, 30);

	mpf_set_d(a, 2.0);
	mpf_set_d(b, 1.0);

	mpf_add(c,a,b);

	gmp_printf("fixed point mof %. *Ff with %d digits\n", 20, c, 20);
	return 0;
}
