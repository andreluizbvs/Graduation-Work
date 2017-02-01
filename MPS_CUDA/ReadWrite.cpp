#include "ReadWrite.h"

ReadWrite::ReadWrite()
{

}

ReadWrite::~ReadWrite()
{

}


ReadWrite::ReadWrite(const ReadWrite &r)
{

}

void ReadWrite::ReadParticles2D(Particle2D* particles,char* name, int nump)
{
	ifstream in;
	in.open(name);

	while(!in.is_open())
	{
		cout<<"trying to read again"<<endl;
		in.open(name);
	}

	in>>nump;

	for(int a=0; a<nump;a++)
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

		particles[a] = particle;
	}
	in.close();
}

void ReadWrite::ReadParticles3D(Particle3D* particles, char* name, int nump)
{
	ifstream in;
	in.open(name);

	while (!in.is_open())
	{
		cout << "trying to read again" << endl;
		in.open(name);
	}

	in >> nump;

	for (int a = 0; a<nump; a++)
	{
		Particle3D particle;

		int type;

		Point3D position, velocity;

		double n, p;

		in >> type;

		in >> position.x;
		in >> position.y;
		in >> position.z;

		in >> velocity.x;
		in >> velocity.y;
		in >> velocity.z;

		in >> p;
		in >> n;

		particle.set_all(position, velocity, type, n, p);

		particles[a] = particle;
	}
	in.close();
}

void ReadWrite::ReadData(char *name, data_in *dados)
{
	ifstream in;

	in.open(name);
	while(!in.is_open())
	{
		cout<<"trying to read again"<<endl;
		in.open(name);
	}

	in>>dados->NumberOfFluid;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->wall_type;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->wall_p_type;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->DensityFluid0;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->DensityFluid1;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->DensityWall_P;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->DensityWall;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->ExpansionFluid0;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->ExpansionFluid1;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->ExpansionWall_P;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->ExpansionWall;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->Gravity;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->MaxTime;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->MaxIteration;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->MaxDt;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->DtRatio;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->Convergence;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->MaxIterPress;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->MinIterPress;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->AveParticleDis;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->rePND;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->reICCG;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->minlen;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->dirichlet;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->col_rat;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in>>dados->lam_rat;
	in.ignore(numeric_limits<streamsize>::max(),'\n');

	in.close();

	dados->set_after_write();
}

void ReadWrite::WriteOut2D(Particle2D* particles, int number, int nump)
{
	FILE*out;
	char nome[1000] = "\0";
	int n = sprintf_s(nome,"saida_tg\\out%d.vtu",number);

	fopen_s(&out,nome,"w");
	while (out == NULL){
		fopen_s(&out, nome, "w");
		cout << "trying to write vtu" << endl;
	}
	fprintf(out, "<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(out, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(out, " <UnstructuredGrid>\n");

	fprintf(out, "  <Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", nump, nump);

	fprintf(out,"   <Points>\n");
	fprintf(out, "    <DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");

	int b;
	for(b = 0; b < nump; b++)
		fprintf(out, "%.15lf %.15lf 0 ", (particles)[b].get_po().x, (particles)[b].get_po().y);

	fprintf(out, "\n");
	fprintf(out,"    </DataArray>\n");
	fprintf(out,"   </Points>\n");

	fprintf(out, "   <PointData>\n");

	fprintf(out,"    <DataArray NumberOfComponents='3' type='Float32' Name='Velocity' format='ascii'>\n");
	for (b = 0; b < nump; b++)
		fprintf(out, "%.15lf %.15lf 0 ", (particles)[b].get_v().x, (particles)[b].get_v().y);


	fprintf(out,"\n");
	fprintf(out, "    </DataArray>\n");

	fprintf(out,"    <DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");


	for (b = 0; b < nump; b++)
		fprintf(out, "%.15lf ", (particles)[b].get_pr());
	
	fprintf(out,"\n");
	fprintf(out,"    </DataArray>\n");

	fprintf(out,"   </PointData>\n");

	//VTK specific information
	fprintf(out,"   <Cells>\n");

	//Connectivity
	fprintf(out,"    <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for (b = 0; b < nump; b++)
	{
		fprintf(out, "%d ",b);
	}
	fprintf(out,"\n");
	fprintf(out,"    </DataArray>\n");

	//Offsets
	fprintf(out,"    <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for (b = 1; b <= nump; b++)
		fprintf(out, "%d ",b);

	fprintf(out,"\n");
	fprintf(out,"    </DataArray>\n");

	//Types
	fprintf(out,"    <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for (b = 0; b < nump; b++)
		fprintf(out,"1 ");

	fprintf(out,"\n");
	fprintf(out,"    </DataArray>\n");

	fprintf(out,"   </Cells>\n");

	fprintf(out,"  </Piece>\n");
	fprintf(out," </UnstructuredGrid>\n");
	fprintf(out,"</VTKFile>");

	fclose(out);
}

void ReadWrite::WriteOut(Particle3D* particles, int number, int nump)
{
	FILE*out;
	char nome[1000] = "\0";
	int n = sprintf_s(nome, "saida_tg\\out%d.vtu", number);

	fopen_s(&out, nome, "w");
	while (out == NULL){
		fopen_s(&out, nome, "w");
		cout << "trying to write vtu" << endl;
	}
	fprintf(out, "<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(out, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(out, " <UnstructuredGrid>\n");

	fprintf(out, "  <Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", nump, nump);

	fprintf(out, "   <Points>\n");
	fprintf(out, "    <DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");

	int b;
	for (b = 0; b < nump; b++)
		fprintf(out, "%.15lf %.15lf %.15lf ", (particles)[b].get_po().x, (particles)[b].get_po().y, (particles)[b].get_po().z);

	fprintf(out, "\n");
	fprintf(out, "    </DataArray>\n");
	fprintf(out, "   </Points>\n");

	fprintf(out, "   <PointData>\n");

	fprintf(out, "    <DataArray NumberOfComponents='3' type='Float32' Name='Velocity' format='ascii'>\n");
	for (b = 0; b < nump; b++)
		fprintf(out, "%.15lf %.15lf %.15lf ", (particles)[b].get_v().x, (particles)[b].get_v().y, (particles)[b].get_po().z);


	fprintf(out, "\n");
	fprintf(out, "    </DataArray>\n");

	fprintf(out, "    <DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");


	for (b = 0; b < nump; b++)
		fprintf(out, "%.15lf ", (particles)[b].get_pr());

	fprintf(out, "\n");
	fprintf(out, "    </DataArray>\n");

	fprintf(out, "   </PointData>\n");

	//VTK specific information
	fprintf(out, "   <Cells>\n");

	//Connectivity
	fprintf(out, "    <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for (b = 0; b < nump; b++)
	{
		fprintf(out, "%d ", b);
	}
	fprintf(out, "\n");
	fprintf(out, "    </DataArray>\n");

	//Offsets
	fprintf(out, "    <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for (b = 1; b <= nump; b++)
		fprintf(out, "%d ", b);

	fprintf(out, "\n");
	fprintf(out, "    </DataArray>\n");

	//Types
	fprintf(out, "    <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for (b = 0; b < nump; b++)
		fprintf(out, "1 ");

	fprintf(out, "\n");
	fprintf(out, "    </DataArray>\n");

	fprintf(out, "   </Cells>\n");

	fprintf(out, "  </Piece>\n");
	fprintf(out, " </UnstructuredGrid>\n");
	fprintf(out, "</VTKFile>");

	fclose(out);
}