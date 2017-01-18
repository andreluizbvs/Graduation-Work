#include "Particles.h"


// estrutura de uma particula 
// uma particula tem posição, velocidade, pressão e dencidade de numero de particulas (dencidade de vizinhança)

Particle2D::Particle2D()
{
	this->position.x=0;
	this->position.y=0;

	this->velocity.x = 0;
	this->velocity.y = 0;
	
	this->Type = 0;
	this->n = 0.0;
	this->p = 0.0;
}

Particle2D::~Particle2D()
{

}

Particle2D::Particle2D(const Particle2D &p)
{
	this->position  = p.position;
	this->velocity = p.velocity;

	this->Type = p.Type;
	this->n = p.n;
	this->p = p.p;
}

void Particle2D::get_all(Point2D &position, Point2D &velocity, int &Type, double &n, double &p)
{
	position = this->position;
	velocity = this->velocity;

	Type = this->Type;
	n = this->n;
	p =this->p;
}

Point2D Particle2D::get_po()
{
	return this->position;
}

Point2D Particle2D::get_v()
{
	return this->velocity;
}

double  Particle2D::get_n()
{
	return this->n;
}

double Particle2D::get_pr()
{
	return this->p;
}

void Particle2D::set_all(Point2D position, Point2D velocity, int Type, double n, double p)
{
	this->position = position;
	this->velocity = velocity;

	this->Type = Type;

	this->n = n;
	this->p = p;
}

void Particle2D::set_n( double n)
{
	this->n = n;
}

void Particle2D::set_po(double x, double y)
{
	this->position.x = x;
	this->position.y = y;
}

void Particle2D::set_pr( double p)
{
	this->p = p;
}

void Particle2D::set_v(double x,double y)
{
	this->velocity.x = x;
	this->velocity.y = y;

}

bool Particle2D::is_border()
{
	if(Type==2)
		return true;
	else
		return false;
}

bool Particle2D::is_fluid()
{
	if(Type==0)
		return true;
	else
		return false;
}


bool Particle2D::is_wall()
{
	if(Type==3)
		return true;
	else
		return false;
}

int Particle2D::get_type()
{
	return this->Type;
}
