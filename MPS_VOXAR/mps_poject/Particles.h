#ifndef __PARTICLES_H__
#define __PARTICLES_H__

using namespace std;

struct Point2D
{
	double x;
	double y;
};

struct Point_3D
{
	double x;
	double y;
	double z;
};

class Particle2D
{
protected:
	Point2D position;
	Point2D velocity;
	
	int Type;
	double n;
	double p;

public:
	Particle2D();//base
	~Particle2D();
	Particle2D(const Particle2D &p);

	void get_all(Point2D &position, Point2D &velocity, int &Type, double &n, double &p);
	Point2D get_po();
	Point2D get_v();
	double get_n();
	double get_pr();
	int get_type();

	void set_all(Point2D position, Point2D velocity, int Type, double n, double p);
	void set_po(double x,double y);
	void set_v(double x,double y);
	void set_n(double n);
	void set_pr(double p);
	
	bool is_fluid();
	bool is_wall();
	bool is_border();
};


#endif