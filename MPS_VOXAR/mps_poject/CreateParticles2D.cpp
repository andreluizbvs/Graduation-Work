#include "Functions.h"
#include "Classes.h"

void createParticles2D(vector<Particle2D>* particles,char* name)
{
	ifstream in;
	in.open(name);

	//FILE*test;
	//fopen_s(&test,"hue.txt","w");
	//fclose(test);

	if(in.is_open())
	{
		cout<<"nop"<<endl;
		system("pause");
		exit(1);
	}

	int num_Part=0;

	in>>num_Part;

	for(int a=0; a<num_Part;a++)
	{
		Particle2D particle;

		int type;

		Point2D position,velocity;

		double n,p;

		in>>type;
		
		in>>position.x;
		in>>position.y;

		in>>velocity.x;
		in>>velocity.y;

		in>>p;
		in>>n;

		particle.set_all(position,velocity,type,n,p);

		particles->push_back(particle);
	}
}
