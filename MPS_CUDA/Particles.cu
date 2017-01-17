#include "Particles.h"

// Estrutura de uma partícula 
// Uma particula tem posição, velocidade, pressão e densidade de número de partículas (densidade de vizinhança)

__host__ __device__ Particle2D::Particle2D()
{
	this->position.x = 0;
	this->position.y = 0;

	this->velocity.x = 0;
	this->velocity.y = 0;

	this->Type = 0;
	this->n = 0.0;
	this->p = 0.0;
}

__host__ __device__ Particle2D::~Particle2D()
{

}

__host__ __device__ Particle2D::Particle2D(const Particle2D &p)
{
	this->position = p.position;
	this->velocity = p.velocity;

	this->Type = p.Type;
	this->n = p.n;
	this->p = p.p;
}

__host__ __device__ void Particle2D::get_all(Point2D &position, Point2D &velocity, int &Type, double &n, double &p)
{
	position = this->position;
	velocity = this->velocity;

	Type = this->Type;
	n = this->n;
	p = this->p;
}

__host__ __device__ Point2D Particle2D::get_po()
{
	return this->position;
}

__host__ __device__ Point2D Particle2D::get_v()
{
	return this->velocity;
}

__host__ __device__ double  Particle2D::get_n()
{
	return this->n;
}

__host__ __device__ double Particle2D::get_pr()
{
	return this->p;
}

__host__ __device__ void Particle2D::set_all(Point2D position, Point2D velocity, int Type, double n, double p)
{
	this->position = position;
	this->velocity = velocity;

	this->Type = Type;

	this->n = n;
	this->p = p;
}

__host__ __device__ void Particle2D::set_n(double n)
{
	this->n = n;
}

__host__ __device__ void Particle2D::set_po(double x, double y)
{
	this->position.x = x;
	this->position.y = y;
}

__host__ __device__ void Particle2D::set_pr(double p)
{
	this->p = p;
}

__host__ __device__ void Particle2D::set_v(double x, double y)
{
	this->velocity.x = x;
	this->velocity.y = y;

}

__host__ __device__ bool Particle2D::is_border()
{
	if (Type == 2)
		return true;
	else
		return false;
}

__host__ __device__ bool Particle2D::is_fluid()
{
	if (Type == 0)
		return true;
	else
		return false;
}

__host__ __device__ bool Particle2D::is_wall()
{
	if (Type == 3)
		return true;
	else
		return false;
}

__host__ __device__ int Particle2D::get_type()
{
	return this->Type;
}


////////////////////////////////////////////////////////////////////////////////////////// 3 DIMENSÕES


__host__ __device__ Particle3D::Particle3D()
{
	this->position.x = 0;
	this->position.y = 0;
	this->position.z = 0;

	this->velocity.x = 0;
	this->velocity.y = 0;
	this->velocity.z = 0;

	this->Type = 0;
	this->n = 0.0;
	this->p = 0.0;
}

__host__ __device__ Particle3D::~Particle3D()
{

}

__host__ __device__ Particle3D::Particle3D(const Particle3D &p)
{
	this->position = p.position;
	this->velocity = p.velocity;

	this->Type = p.Type;
	this->n = p.n;
	this->p = p.p;
}

__host__ __device__ void Particle3D::get_all(Point3D &position, Point3D &velocity, int &Type, double &n, double &p)
{
	position = this->position;
	velocity = this->velocity;

	Type = this->Type;
	n = this->n;
	p = this->p;
}

__host__ __device__ Point3D Particle3D::get_po()
{
	return this->position;
}

__host__ __device__ Point3D Particle3D::get_v()
{
	return this->velocity;
}

__host__ __device__ double  Particle3D::get_n()
{
	return this->n;
}

__host__ __device__ double Particle3D::get_pr()
{
	return this->p;
}

__host__ __device__ void Particle3D::set_all(Point3D position, Point3D velocity, int Type, double n, double p)
{
	this->position = position;
	this->velocity = velocity;

	this->Type = Type;

	this->n = n;
	this->p = p;
}

__host__ __device__ void Particle3D::set_n(double n)
{
	this->n = n;
}

__host__ __device__ void Particle3D::set_po(double x, double y, double z)
{
	this->position.x = x;
	this->position.y = y;
	this->position.z = z;
}

__host__ __device__ void Particle3D::set_pr(double p)
{
	this->p = p;
}

__host__ __device__ void Particle3D::set_v(double x, double y, double z)
{
	this->velocity.x = x;
	this->velocity.y = y;
	this->velocity.z = z;

}

__host__ __device__ bool Particle3D::is_border()
{
	if (Type == 2)
		return true;
	else
		return false;
}

__host__ __device__ bool Particle3D::is_fluid()
{
	if (Type == 0)
		return true;
	else
		return false;
}

__host__ __device__ bool Particle3D::is_wall()
{
	if (Type == 3)
		return true;
	else
		return false;
}

__host__ __device__ int Particle3D::get_type()
{
	return this->Type;
}