#include "Kernel.h"

double kernel::weight(double r,double re)
{
	double ORIG_rr;
	double ORIG_ww;

	ORIG_rr = (r / re);
	if (ORIG_rr<1.0) ORIG_ww = 1.0 / ORIG_rr - 1.0;
	else ORIG_ww = 0.0;

	return(ORIG_ww);
}

kernel::kernel()
{

}

kernel::~kernel()
{

}

kernel::kernel(const kernel &k)
{

}