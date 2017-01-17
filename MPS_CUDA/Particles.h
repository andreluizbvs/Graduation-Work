#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

using namespace std;

template <typename T>
struct Vector2D
{
	T x;
	T y;
};

struct Point2D
{
	double x;
	double y;
	//double z;
};

struct Point3D
{
	double x;
	double y;
	double z;
};

class Particle2D
{
public:
//protected:
	Point2D position;
	Point2D velocity;
	
	int Type;
	double n;
	double p;

//public:
	__host__ __device__ Particle2D();//base
	__host__ __device__ ~Particle2D();
	__host__ __device__ Particle2D(const Particle2D &p);

	__host__ __device__ void get_all(Point2D &position, Point2D &velocity, int &Type, double &n, double &p);
	__host__ __device__ Point2D get_po();
	__host__ __device__ Point2D get_v();
	__host__ __device__ double get_n();
	__host__ __device__ double get_pr();
	__host__ __device__ int get_type();

	__host__ __device__ void set_all(Point2D position, Point2D velocity, int Type, double n, double p);
	__host__ __device__ void set_po(double x, double y);
	__host__ __device__ void set_v(double x, double y);
	__host__ __device__ void set_n(double n);
	__host__ __device__ void set_pr(double p);
	
	__host__ __device__ bool is_fluid();
	__host__ __device__ bool is_wall();
	__host__ __device__ bool is_border();
};

class Particle3D
{
public:
	//protected:
	Point3D position;
	Point3D velocity;

	int Type;
	double n;
	double p;

	//public:
	__host__ __device__ Particle3D();//base
	__host__ __device__ ~Particle3D();
	__host__ __device__ Particle3D(const Particle3D &p);

	__host__ __device__ void get_all(Point3D &position, Point3D &velocity, int &Type, double &n, double &p);
	__host__ __device__ Point3D get_po();
	__host__ __device__ Point3D get_v();
	__host__ __device__ double get_n();
	__host__ __device__ double get_pr();
	__host__ __device__ int get_type();

	__host__ __device__ void set_all(Point3D position, Point3D velocity, int Type, double n, double p);
	__host__ __device__ void set_po(double x, double y, double z);
	__host__ __device__ void set_v(double x, double y, double z);
	__host__ __device__ void set_n(double n);
	__host__ __device__ void set_pr(double p);

	__host__ __device__ bool is_fluid();
	__host__ __device__ bool is_wall();
	__host__ __device__ bool is_border();
};

#endif