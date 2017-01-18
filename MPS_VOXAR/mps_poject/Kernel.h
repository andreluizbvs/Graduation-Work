#ifndef __KERNEL_H__
#define __KERNEL_H__

class kernel
{
public:
	double weight(double r,double re);
	
	kernel();
	~kernel();
	kernel(const kernel &k);

};

#endif