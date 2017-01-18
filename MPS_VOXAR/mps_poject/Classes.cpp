#include "Classes.h"

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

void Particle2D::get_all(Point2D &position, Point2D &velocity, int &Type, double &n, double &p)
{
	position = this->position;
	velocity = this->velocity;

	Type = this->Type;
	n = this->n;
	p =this->p;
}

void Particle2D::set_all(Point2D position, Point2D velocity, int Type, double n, double p)
{
	this->position = position;
	this->velocity = velocity;

	this->Type = Type;

	this->n = n;
	this->p = p;
}

