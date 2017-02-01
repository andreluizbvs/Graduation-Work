/*

    MPS Code "mps.c"

*/
/*	MPS : Moving Particle Semi-implicit
		This is a particle method which can be applied to
		incompressible flow without grids. Particle motion
		is just fluid motion; this is fully Lagragian calculation.
		The MPS method was invented by S. Koshizuka, University
		of Tokyo. If anyone who has questions about this code
		can contact him at koshizuka@q.t.u-tokyo.ac.jp .
*/
/*	mps.c : A fundamental 2-D code of the MPS method.
			Viscosity and surface tension terms are not included.
*/
/*	#### DO NOT USE THIS CODE FOR COMMERCIAL ACTIVITIES. #### */

/*  copyright (c) 2005 by Seiichi Koshizuka */

/* #include <cstdlib> */
/* This line is necessary for C++ compiler.
                      Remove this line in C compiler */
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <string.h>

#define YES 0
#define NO 1

#define NN 3000
#define DIM 2
#define FLUID 5
#define IN_PARAM "mps.data"
#define IN_PROF "mps.grid"
#define OUT_PROF "mps.prof"
#define saida_tela "tela.out"
#define NEIGHBOR 60
#define VLIM 0.1
#define DIAG 2.0
#define GHOST -1
#define DT_TOO_SMALL 0.05
#define INCREASE_DT 1.1
#define FILE_OUTPUT_LIMIT 10000000000
#define MOD_LIMIT 3

void cal_n();
double weight();
double cal_dt();
void cal_buoyancy();
void cal_convection();
void set_source();
void set_matrix();
void set_bcon();
void check_bcon();
void incom_decomp();
int pcg_solver();
void solver_ll();
void rev_pressgrad();
void rev_convection();
void reset_neigh();
void set_neigh();
void file1input();
void file2input();
void fileoutput();
void timeovercheck();
FILE *filecheckopen();
void set_vector_zero();
void copy_vectors();
void add_vectors();
void sub_vectors();
void mul_matrix_vector();
void mul_vectors();
void collision();
int checkconv();
int outputcheck();
int outputcount=0;
void exit_job();
void set_pmin();
int dt_too_small();
double rho2pg();
double rho2sm();
FILE*tela;
int nump=-1; /* number of all particles */
int itmx=-1; /* maximum time iteration */
char nome[10000] = {'\0'};

typedef struct Ponto
{
    double x;
    double y;
} Ponto;

void converte_saida()
{
    printf("%d",itmx);
    FILE* entrada = fopen( OUT_PROF,"r" );
    double time_step;

    FILE *output;
    Ponto *locais;
    double **vx, *presao;
    int lixo;
    double N;
    int i;
    int ntotal = nump;

    for(i=1; i<=itmx; i++)
    {
        sprintf (nome, "saida_navegante//out%d.vtu", i);
        FILE* saida_navegante = fopen( nome,"w" );

        fscanf(entrada," %lf",&time_step);
        fscanf(entrada," %d",&lixo);

        locais = (Ponto*) malloc(sizeof(Ponto) * nump);
        vx = (double**) malloc(sizeof(double*) * 2);
        vx[0] = (double*) malloc(sizeof(double) * nump);
        vx[1] = (double*) malloc(sizeof(double) * nump);
        presao = (double*) malloc(sizeof(double) * nump);


        int a;
        for(a=0; a<nump; a++)
        {
            fscanf(entrada," %d",&lixo);
            fscanf(entrada," %lf",&(locais[a].x));
            fscanf(entrada," %lf",&(locais[a].y));
            fscanf(entrada," %lf",&(vx[0][a]));
            fscanf(entrada," %lf",&(vx[1][a]));
            fscanf(entrada," %lf",&presao[a]);
            fscanf(entrada," %lf",&N);
        }

        /*******************************************************************/
        fprintf(saida_navegante, "<?xml version='1.0' encoding='UTF-8'?>\n");
        fprintf(saida_navegante, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
        fprintf(saida_navegante, " <UnstructuredGrid>\n");

        fprintf(saida_navegante,"  <Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",ntotal,ntotal);

        fprintf(saida_navegante,"   <Points>\n");
        fprintf(saida_navegante, "    <DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");

        int b;
        for(b = 0; b < ntotal; b++)
            fprintf(saida_navegante, "%lf %lf 0 ", locais[b].x, locais[b].y);

        fprintf(saida_navegante, "\n");
        fprintf(saida_navegante,"    </DataArray>\n");
        fprintf(saida_navegante,"   </Points>\n");

        fprintf(saida_navegante, "   <PointData>\n");



        fprintf(saida_navegante,"    <DataArray NumberOfComponents='3' type='Float32' Name='Velocity' format='ascii'>\n");
        for(b = 0; b < ntotal; b++)
            fprintf(saida_navegante, "%lf %lf 0 ", vx[0][b], vx[1][b]);


        fprintf(saida_navegante,"\n");
        fprintf(saida_navegante, "    </DataArray>\n");

        fprintf(saida_navegante,"    <DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");


        for( b=0; b < ntotal; b++)
        {
            fprintf(saida_navegante, "%lf ", presao[b]);
        }
        fprintf(saida_navegante,"\n");
        fprintf(saida_navegante,"    </DataArray>\n");

        //More to add soon...

        fprintf(saida_navegante,"   </PointData>\n");

        //VTK specific information
        fprintf(saida_navegante,"   <Cells>\n");

        //Connectivity
        fprintf(saida_navegante,"    <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
        for( b = 0; b < ntotal; b++)
        {
            fprintf(saida_navegante, "%d ",b);
        }
        fprintf(saida_navegante,"\n");
        fprintf(saida_navegante,"    </DataArray>\n");

        //Offsets
        fprintf(saida_navegante,"    <DataArray type='Int32' Name='offsets' format='ascii'>\n");
        for( b = 1; b <= ntotal; b++)
            fprintf(saida_navegante, "%d ",b);

        fprintf(saida_navegante,"\n");
        fprintf(saida_navegante,"    </DataArray>\n");

        //Types
        fprintf(saida_navegante,"    <DataArray type='UInt8' Name='types' format='ascii'>\n");
        for(b = 0; b < ntotal; b++)
            fprintf(saida_navegante,"1 ");

        fprintf(saida_navegante,"\n");
        fprintf(saida_navegante,"    </DataArray>\n");

        fprintf(saida_navegante,"   </Cells>\n");

        fprintf(saida_navegante,"  </Piece>\n");
        fprintf(saida_navegante," </UnstructuredGrid>\n");
        fprintf(saida_navegante,"</VTKFile>");



        fclose(saida_navegante);

        free(locais);
        free(vx[0]);
        free(vx[1]);
        free(vx);
        free(presao);
    }
}

int main()
{
    int fluid=-1; /* number of fluid type */
    double rho[FLUID]= {-1}; /* density */
    double alpha[FLUID]= {-1}; /* compressibility */
    double grav=-1; /* gravitational constant */
    double tmax=-1,dtmx=-1,dtra=-1,epsp=-1; /* computation parameters */
    double disa=-1,rerp=-1,rericcg=-1; /* kernel parameters */
    int imxp1=-1,imnp1=-1; /* computation parameters */
    char outputtype[3]; /* output time; I : iteration, T : time */
    int interval=-1; /* output interval */
    double output_dt=-1, output_time=-1; /* output interval, output time */
    double time_sim=-1; /* simulation time [sec] */
    static double x[DIM][NN],v[DIM][NN],p[NN];
    static double xn[DIM][NN],vn[DIM][NN],dv[DIM][NN],dp[NN],source[NN];
    static double n[NN],n_fluid[FLUID][NN];
    static int typep[NN]; /* particle type */
    double rep=-1,rep2=-1,reiccg=-1,reiccg2=-1; /* re of the kernel */
    double n0p=-1,n0iccg=-1; /* fixed particle number density */
    int n0p_particle=-1; /* particle to evaluate n0p in the initial configuration */
    double dt=-1;
    double lambda=-1; /* int (w r^2 dv) / int (w dv) */
    int it=-1;
    static int neigh[NN][NEIGHBOR]; /* neighboring particle numbers storage */
    static int neigh_iccg[NN][NEIGHBOR];
    static double poiss[NN][NEIGHBOR],ic[NN][NEIGHBOR];
    static int bcon[NN]; /* -1 : ignore (wall,symmetry), 0 : to be solved, 1 : value fixed */
    static int check[NN];
    double r2lim=-1;
    int i=-1,j=-1,k=-1,l=-1,kk=-1;
    FILE *fp1,*fp2,*fpo;
    double minlen=-1,dirichlet=-1,col_rat=-1;
    int num_dirichlet=0;
    double ra=-1;
    int wall_type=-1,wall_p_type=-1;
    static double pmin[NN];
    double lam_rat=-1;
    tela = fopen("saida.out","w");

    /*                */
    /*  File setting  */
    /*                */

    fprintf(tela,"%%%% mps %%%%\n");

    /*  filename setting and open  % a.out fname1 fname2 fnameo */

    fp1=filecheckopen(1,IN_PARAM,"r");
    fp2=filecheckopen(2,IN_PROF,"r");
    fpo=filecheckopen(3,OUT_PROF,"w");

    /*  data input from fp1 */
    /*  input from %%%%.data */

    file1input(fp1,&fluid,&wall_type,&wall_p_type,rho,alpha,&grav,
               &tmax,&itmx,&dtmx,&dtra,&epsp,&imxp1,&imnp1,outputtype,&interval,&output_dt,&disa,&rerp,&rericcg,
               &minlen,&dirichlet,&col_rat,&lam_rat,
               &n0p,&n0iccg,&n0p_particle);

    /*  data input from fp2, initial profile */
    /*  input from %%%%.grid */

    file2input(fp2,&time_sim,&nump,&typep,x,v,p,n);

    /*                                 */
    /*  calculation parameters setting */
    /*                                 */
    rep=disa*rerp; /* real distance [m] of re(PND) */
    rep2=rep*rep;
    reiccg=disa*rericcg; /* real distance [m] of re(ICCG) */
    reiccg2=reiccg*reiccg;

    /* If n0p=0.0 in the input data, n0p is given as n[n0p_particle] */
    /* "n0p_particle" is given in the input data */

    if(n0p==0.0)
    {
        reset_neigh(nump,neigh);
        reset_neigh(nump,neigh_iccg);
        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST ) continue;
            set_neigh(i,x,nump,typep,neigh,rep2,neigh_iccg,reiccg2,wall_type);
        }
        cal_n(n0p_particle,x,n,n_fluid,neigh,typep,rep,fluid,wall_type,wall_p_type);
        n0p=n[n0p_particle];
        cal_n(n0p_particle,x,n,n_fluid,neigh_iccg,typep,reiccg,fluid,wall_type,wall_p_type);
        n0iccg=n[n0p_particle];
        fprintf(tela,"calculated n0p, n0iccg = %e, %e\n",n0p,n0iccg);
    }

    /* Particles collide when r < r2lim */

    r2lim=disa*disa*minlen*minlen;

    /* pi * ra^2 = disa^2 */

    ra=disa/1.7724539;

    lambda=lam_rat*(reiccg2*reiccg2/12.0-reiccg*ra*ra*ra/3.0+ra*ra*ra*ra/4.0)
           /(reiccg2/2.0-reiccg*ra+ra*ra/2.0);

    fprintf(tela,"rep,reiccg=%e,%e\n",rep,reiccg);
    fprintf(tela,"r2lim=%e\n",r2lim);
    fprintf(tela,"lambda=%e\n",lambda);

    /* Initla dt setting */

    dt=dtmx;

    /* Initial setting of neighboring tables */

    reset_neigh(nump,neigh);
    reset_neigh(nump,neigh_iccg);
    for(i=0; i<nump; i++)
    {
        if ( typep[i]==GHOST ) continue;
        set_neigh(i,x,nump,typep,neigh,rep2,neigh_iccg,reiccg2,wall_type);
        cal_n(i,x,n,n_fluid,neigh,typep,rep,fluid,wall_type,wall_p_type);
    }

    /* Initial setting of variables */

    copy_vectors(nump,x[0],xn[0]);
    copy_vectors(nump,x[1],xn[1]);
    copy_vectors(nump,v[0],vn[0]);
    copy_vectors(nump,v[1],vn[1]);

    /* Next output timing */

    if(outputtype[0]=='T')output_time=time_sim+output_dt;

    /* Initial particles' profile is output to %%%%.prof file */

    fileoutput(fpo,time_sim,nump,typep,xn,vn,p,n);

    /*             */
    /*  MAIN LOOP  */
    /*             */

    for(it=1; it<=itmx; it++)
    {

        /*  dt control and time advancing */

        dt=cal_dt(typep,v,nump,dtra,disa,it,dtmx,wall_type,wall_p_type,dt,INCREASE_DT);

        time_sim=time_sim+dt;

        //fprintf(tela,">> it,dt,time_sim,type=%d,%f,%f,%d\n",it,dt,time_sim,i);

        /* If dt is too small, the run stops */

        if( dt_too_small(dt,dtmx,DT_TOO_SMALL)==YES )
        {
            fprintf(tela,"**** The time step falls to a too small value, calculaion ends.  dt = %f\n",dt);
            //exit_job();
        }

        /*  buoyancy force */

        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST || typep[i]==wall_type || typep[i]==wall_p_type ) continue;
            cal_buoyancy(i,vn,grav,dt);
        }

        /*  convection (particle motion) */

        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST || typep[i]==wall_type || typep[i]==wall_p_type ) continue;
            cal_convection(i,xn,vn,dt);
        }

        /*  collision */
        /* This is for avoiding particles' penetration into walls */
        /* If this function works, message is output */
        /* r2lim is the main parameter */

        collision(nump,typep,neigh,xn,vn,rho,r2lim,dt,rep2,col_rat,wall_type,wall_p_type);

        /*  setting of neighboring tables */

        reset_neigh(nump,neigh);
        reset_neigh(nump,neigh_iccg);
        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST ) continue;
            set_neigh(i,xn,nump,typep,neigh,rep2,neigh_iccg,reiccg2,wall_type);
            /*    fprintf(tela,"i,neigh=%d,%d\n",i,neigh[i][neigh[i][0]]); */
        }

        /*                       */
        /*  pressure correction  */
        /*        ICCG           */

        /*  particle number density */

        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST || typep[i]==wall_type ) continue;
            cal_n(i,xn,n,n_fluid,neigh,typep,rep,fluid,wall_type,wall_p_type);
            /*        fprintf(tela,"i,n=%d,%e\n",i,n[i]); */
        }

        /*  boundary codition setting */

        for(i=0; i<nump; i++)
        {
            set_bcon(i,typep,bcon,n,n0p,dirichlet,wall_type);
            /*            fprintf(tela,"i,bcon=%d,%d\n",i,bcon[i]);  */
        }

        /*  source term */

        num_dirichlet=0;
        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST || typep[i]==wall_type ) continue;
            if(bcon[i]==0)set_source(i,typep,source,n,n0p,rho,dt);
            else if(bcon[i]==1)
            {
                source[i]=0.0;
                if(typep[i]==0)num_dirichlet++;
            }
            /*               fprintf(tela,"i,source=%d,%e\n",i,source[i]); */
        }

        check_bcon(nump,typep,bcon,neigh,check);
        /*******************************************************************************************************/
        /*  coefficient matrix */

        for(i=0; i<nump; i++)
        {
            set_matrix(i,typep,bcon,neigh_iccg,poiss,xn,n0iccg,reiccg,reiccg2,lambda,rho,alpha,dt);
        }

        for(i=0; i<nump; i++)
        {
            //fprintf(tela,"i,poiss=%d,%e\n",i,poiss[i][0]);
        }



        /*  coefficient check */

        for(i=0; i<nump; i++)
        {
            if(check[i]==0)
            {
                //fprintf(tela,"i,typep,bcon,check,n,source,poiss,neigh = %d,%d,%d,%d,%f,%f,%f,%d\n",i,typep[i],bcon[i],check[i],n[i],source[i],poiss[i][0],neigh[i][0]);
                for(l=1; l<= neigh[i][0]; l++)
                {
                    j=neigh[i][l];
                    //fprintf(tela,"     l,j,typep,bcon,check,n,poiss = %d,%d,%d,%d,%d,%f,%f\n",l,j,typep[i],bcon[j],check[j],n[j],poiss[i][l]);
                    for(kk=1; kk<=neigh[j][0]; kk++)
                    {
                        k=neigh[j][kk];
                        fprintf(tela,"%d,",k);
                    }
                    fprintf(tela,"\n");
                }
                poiss[i][0]=poiss[i][0]*DIAG;
            }
        }

        /*  incomplete Choleskey decomposition */

        incom_decomp(nump,bcon,neigh_iccg,poiss,ic);

        for(i=0; i<nump; i++)
        {
            /*fprintf(tela,"i,ic=%d,%e\n",i,ic[i][0])*/;
        }

        /*  preconditioned conjugate gradient solver */

        set_vector_zero(nump,dp);

        /* ICCG iteration */

        pcg_solver(nump,bcon,neigh_iccg,poiss,dp,source,ic,epsp,imnp1,imxp1,dt);

        /*  minus pressure -> zero pressure */

        for(i=0; i<nump; i++)
        {
            if(dp[i]<0.0)dp[i]=0.0;
            //fprintf(tela,"i,dp=%d,%e\n",i,dp[i]);
        }

        /*  maximum pressure table generation */

        for(i=0; i<nump; i++)
        {
            if ( bcon[i]==-1 || typep[i]==wall_p_type ) continue;
            set_pmin(i,typep,neigh_iccg,dp,pmin,wall_type);
            //fprintf(tela,"i,pmin=%d,%e\n",i,pmin[i]);
        }

        /*  velocity correction, pressure gradient term */

        set_vector_zero(nump,dv[0]);
        set_vector_zero(nump,dv[1]);
        for(i=0; i<nump; i++)
        {
            if ( bcon[i]==-1 || typep[i]==wall_p_type ) continue;
            rev_pressgrad(i,typep,neigh_iccg,n,rho,xn,dv,dp,dt,rep2,rep,n0p,wall_type,wall_p_type,pmin);
            //fprintf(tela,"i,dv=%d,%e,%e\n",i,dv[0][i],dv[1][i]);
        }

        add_vectors(nump,vn[0],vn[0],dv[0],1.0);
        add_vectors(nump,vn[1],vn[1],dv[1],1.0);
        for(i=0; i<nump; i++)
        {
            //fprintf(tela,"i,vn=%d,%e,%e\n",i,vn[0][i],vn[1][i]);
        }
        copy_vectors(nump,dp,p);

        /*  correction of particle motion */

        for(i=0; i<nump; i++)
        {
            if ( bcon[i]==-1 || typep[i]==wall_p_type ) continue;
            cal_convection(i,xn,dv,dt);
            //fprintf(tela,"i,xn=%d,%e,%e\n",i,xn[0][i],xn[1][i]);
        }

        /*  output */

        if( outputcheck(outputtype,time_sim,&output_time,output_dt,it,interval) == YES )
            fprintf(tela,"outputcheck =%d,%d\n",outputcheck(outputtype,time_sim,&output_time,output_dt,it,interval),it);
        fileoutput(fpo,time_sim,nump,typep,xn,vn,p,n);

        /*  setting of neighboring tables */

        reset_neigh(nump,neigh);
        reset_neigh(nump,neigh_iccg);
        for(i=0; i<nump; i++)
        {
            if ( typep[i]==GHOST ) continue;
            set_neigh(i,xn,nump,typep,neigh,rep2,neigh_iccg,reiccg2,wall_type);
        }

        /*  copy xn -> x */
        copy_vectors(nump,xn[0],x[0]);
        copy_vectors(nump,xn[1],x[1]);
        copy_vectors(nump,vn[0],v[0]);
        copy_vectors(nump,vn[1],v[1]);

        timeovercheck(time_sim,tmax);

    } /* end of time loop */

    fclose(fpo);
    fclose(tela);

    converte_saida();

    return 0;
}

/*************************************************************************/
/*                                                                       */
/*                     functions                                         */
/*                                                                       */
/*************************************************************************/

/* particle number density */
void
cal_n(i,x,n,n_fluid,neigh,typep,rep,fluid,wall_type,wall_p_type)
int i;
double x[DIM][NN];
double n[NN],n_fluid[FLUID][NN];
int neigh[NN][NEIGHBOR];
int typep[NN];
double rep;
int fluid,wall_type,wall_p_type;
{
    int j,kj,l;
    double xx,yy,zz,nn;

    n[i]=0.0;
    for(j=0; j<fluid; j++)n_fluid[j][i]=0.0;

    for(l=1; l<=neigh[i][0]; l++)
    {
        j=neigh[i][l];
        kj=typep[j];

        xx=(x[0][j]-x[0][i]);
        yy=(x[1][j]-x[1][i]);
        zz=sqrt(xx*xx+yy*yy);

        nn=weight(zz,rep);
        n_fluid[kj][i]+=nn;

    }

    for(j=0; j<fluid; j++)n[i]+=n_fluid[j][i];

}

/* weight function */
double
weight(r,re)
double r,re;
{
    double rr;
    double ww;

    rr=r/re;
    if (rr<1.0) ww=1.0/rr-1.0;
    else ww=0.0;

    return(ww);
}

/* dt automatic setting */
double
cal_dt(typep,v,nump,dtra,disa,it,dtmx,wall_type,wall_p_type,dt_old,ratio)
int typep[NN];
double v[DIM][NN];
int nump;
double dtra,disa;
int it;
double dtmx;
int wall_type,wall_p_type;
double dt_old,ratio;
{
    double dt,dt_lim;
    double vv,vvmax;
    int i;

    vvmax=0.0;
    for(i=0; i<nump; i++)
    {
        if(typep[i]==GHOST || typep[i]==wall_type || typep[i]==wall_p_type)continue;
        vv=v[0][i]*v[0][i]+v[1][i]*v[1][i];
        if(vv>vvmax)vvmax=vv;
    }
    vvmax=sqrt(vvmax);

    dt_lim=dt_old*ratio;
    if(vvmax==0.0)return(dt_lim);

    dt=dtra*disa/vvmax;

    if( dt>dt_lim ) dt=dt_lim;
    if( dt>dtmx ) dt=dtmx;

    return(dt);
}

/* buoyancy force calculation */
void
cal_buoyancy(i,v,grav,dt)
int i;
double v[DIM][NN];
double grav,dt;
{
    v[1][i]=v[1][i]-grav*dt;
}

/* convection, particle motion */
void
cal_convection(i,x,v,dt)
int i;
double x[DIM][NN],v[DIM][NN];
double dt;
{
    x[0][i]=x[0][i]+v[0][i]*dt;
    x[1][i]=x[1][i]+v[1][i]*dt;
}

/* resetting a neighboring table */
void
reset_neigh(nump,neigh)
int nump;
int neigh[NN][NEIGHBOR];
{
    int i;

    for(i=0; i<nump; i++)neigh[i][0]=0;
}

/* setting a neighboring table */
void
set_neigh(i,x,nump,typep,neigh,rep2,neigh_iccg,reiccg2,wall_type)
int i;
double x[DIM][NN];
int nump;
int typep[NN];
int neigh[NN][NEIGHBOR];
double rep2;
int neigh_iccg[NN][NEIGHBOR];
double reiccg2;
int wall_type;
{
    int j;
    double xx,yy;

    for(j=i+1; j<nump; j++)
    {
        if( typep[j]==GHOST )continue;
        if( typep[i]==wall_type && typep[j]==wall_type )continue;

        xx=(x[0][j]-x[0][i]);
        xx=xx*xx;
        if(xx>reiccg2)continue;
        yy=(x[1][j]-x[1][i]);
        xx+=yy*yy;
        if(xx>reiccg2)continue;
        neigh_iccg[i][0]++;
        neigh_iccg[j][0]++;
        if(neigh_iccg[i][0]>=NEIGHBOR)
        {
            fprintf(tela,"#### neighboring particle number above the limit\n");
            fprintf(tela,"i=%d,neigh_iccg[i][0]=%d\n",i,neigh_iccg[i][0]);
            for(j=1; j<=neigh_iccg[i][0]; j++)fprintf(tela," %d ",neigh_iccg[i][j]);
            fprintf(tela," end\n");
            //exit_job();
        }
        neigh_iccg[i][neigh_iccg[i][0]]=j;
        neigh_iccg[j][neigh_iccg[j][0]]=i;

        if(xx>rep2)continue;
        neigh[i][0]++;
        neigh[j][0]++;
        neigh[i][neigh[i][0]]=j;
        neigh[j][neigh[j][0]]=i;

    }

}

FILE
*filecheckopen(i,fname,st)
int i;
char *fname;
char st[];
{
    FILE *fp;

    if ((fp=fopen(fname,st))==NULL)
    {
        fprintf(tela,"file open %s not succeeded\n",fname);
        exit_job();
    }

    return(fp);
}

void
file1input(fp1,fluid,wall_type,wall_p_type,rho,alpha,grav,
           tmax,itmx,dtmx,dtra,epsp,imxp1,imnp1,outputtype,interval,output_dt,disa,rerp,rericcg,
           minlen,dirichlet,col_rat,lam_rat,
           n0p,n0iccg,n0p_particle)
FILE *fp1;
int *fluid;
int *wall_type,*wall_p_type;
double rho[FLUID];
double alpha[FLUID];
double *grav;
double *tmax;
int *itmx;
double *dtmx;
double *dtra;
double *epsp;
int *imxp1,*imnp1;
char outputtype[3];
int *interval;
double *output_dt;
double *disa;
double *rerp,*rericcg;
double *minlen,*dirichlet,*col_rat,*lam_rat;
double *n0p,*n0iccg;
int *n0p_particle;
{
    int i;
    char nomes[100];

    /* physical properties */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %d",fluid); /* number of fluids */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %d",wall_type); /* type of fluid treated as wall */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %d",wall_p_type); /* type of fluid treated as wall+pressure */
    for(i=0; i<*fluid; i++)
    {
        fscanf(fp1,"%s",nomes);
        fscanf(fp1,"%*s %lf",&rho[i]); /* density of water [kg/m3] */
    }
    for(i=0; i<*fluid; i++)
    {
        fscanf(fp1,"%s",nomes);
        fscanf(fp1,"%*s %lf",&alpha[i]); /* expansion coefficient */
    }
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",grav); /* gravitational constant [m/s2] */

    fprintf(tela,"fluid=%d\n",*fluid);
    fprintf(tela,"wall_type,wall_p_type=%d,%d\n",*wall_type,*wall_p_type);
    fprintf(tela,"rho=");
    for(i=0; i<*fluid; i++)fprintf(tela,"%lf ",rho[i]);
    fprintf(tela,"\n");
    fprintf(tela,"alpha=");
    for(i=0; i<*fluid; i++)fprintf(tela,"%e ",alpha[i]);
    fprintf(tela,"\n");
    fprintf(tela,"grav=%lf\n",*grav);

    /*  calculation parameters */

    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",tmax); /* maximum time [sec] */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %d",itmx); /* maximum time iteration */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",dtmx); /* maximum dt [sec] */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",dtra); /* dt ratio to Courant number */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",epsp); /* convergence condition of pressure routine */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %d",imxp1); /* maximum iteration of pressure routine */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %d",imnp1); /* minimum iteration of pressure routine */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %s",outputtype); /* output type; I : iteration, T : time */
    if (outputtype[0]=='I')
    {
        fscanf(fp1,"%s",nomes);
        fscanf(fp1,"%*s %d",interval);
    } /* minimum iteration of pressure routine */
    else
    {
        fscanf(fp1,"%s",nomes);
        fscanf(fp1,"%*s %lf",output_dt); /* output time intererval */
    }
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",disa); /* average particle distance [m] */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",rerp); /* re ratio to disa for continuity equation */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",rericcg); /* re ratio to disa for diffusion terms */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",minlen); /* minimum distance for collision */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",dirichlet); /* dirichlet condition by n/n0 */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",col_rat); /* collision coefficient */
    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",lam_rat); /* relaxation factor in pressure calculation */

    fprintf(tela,"tmax,itmx=%lf,%d\n",*tmax,*itmx);
    fprintf(tela,"dtmx,dtra=%lf,%lf\n",*dtmx,*dtra);
    fprintf(tela,"epsp=%e\n",*epsp);
    fprintf(tela,"imxp1,imnp1,interval=%d,%d\n",*imxp1,*imnp1);
    fprintf(tela,"output : type, i, t = %s,%d,%lf\n",outputtype,*interval,*output_dt);
    fprintf(tela,"disa,rerp,rericcg=%lf,%lf,%lf\n",*disa,*rerp,*rericcg);
    fprintf(tela,"minlen,dirichlet,col_rat,lam_rat=%lf,%lf,%lf,%lf\n",*minlen,*dirichlet,*col_rat,*lam_rat);

    /*  optional parameters */

    fscanf(fp1,"%s",nomes);
    fscanf(fp1,"%*s %lf",n0p); /* fixed PND, If this is 0.0, n(n0p_particle) is used. */
    if(*n0p==0.0)
        fscanf(fp1,"%d",n0p_particle); /* particle number to calculate n0p */
    else
        fscanf(fp1,"%lf",n0iccg); /* fixed PND for iccg */

    fprintf(tela,"n0p=%lf\n",*n0p);
    if(*n0p==0.0)
        fprintf(tela,"n0p_particle=%d\n",*n0p_particle);
    else
        fprintf(tela,"n0iccg=%lf\n",n0iccg);

    fclose(fp1);
}

void
file2input(fp2,time_sim,nump,typep,x,v,p,n)
FILE *fp2;
double *time_sim;
int *nump;
int typep[NN];
double x[DIM][NN],v[DIM][NN],p[NN],n[NN];
{
    int i;

    fscanf(fp2,"%lf",time_sim);
    fprintf(tela,"time_sim=%f\n",*time_sim);
    fscanf(fp2,"%d",nump);
    fprintf(tela,"nump=%d\n",*nump);

    if(*nump>NN)
    {
        fprintf(tela,"#### Number of particles is larger than the limit NN ####\n");
        fprintf(tela,"        nump,NN=%d,%d\n",nump,NN);
        //exit_job();
    }

    for(i=0; i<*nump; i++)
        fscanf(fp2,"%d %lf %lf %lf %lf %lf %lf",&typep[i],&x[0][i],&x[1][i],&v[0][i],&v[1][i],&p[i],&n[i]);

    fclose(fp2);
}

void
timeovercheck(time_sim,tmax)
double time_sim,tmax;
{
    if(time_sim > tmax)
    {
        fprintf(tela,"$$$$ end of calculation, time_sim over, %f\n",time_sim);
        //exit_job();
    }
}

void
fileoutput(fpo,time_sim,nump,typep,x,v,p,n)
FILE *fpo;
double time_sim;
int nump;
int typep[NN];
double x[DIM][NN],v[DIM][NN];
double p[NN],n[NN];
{
    int i;

    fprintf(fpo,"%f\n",time_sim);
    fprintf(fpo,"%d\n",nump);
    for(i=0; i<nump; i++)
        fprintf(fpo,"%d %f %f %f %f %f %f\n",typep[i],x[0][i],x[1][i],v[0][i],v[1][i],p[i],n[i]);//tipo de particula, Pos_X, Pos_y, V_X, V_Y,
    fflush(fpo);

    outputcount++;
    //fprintf(tela,"data output : %d\n",outputcount);

    if(outputcount>=FILE_OUTPUT_LIMIT)
    {
        //fprintf(tela,"data output overflow\n");
        //exit_job();
    }
}

void
copy_vectors(nump,x1,x2)
int nump;
double x1[NN],x2[NN];
{
    int i;

    for(i=0; i<nump; i++) x2[i]=x1[i];
}

void
set_vector_zero(nump,x)
int nump;
double x[NN];
{
    int i;

    for(i=0; i<nump; i++)x[i]=0.0;
}

/* pressure gradient term calculation */
void
rev_pressgrad(i,typep,neigh,n,rho,x,dv,dp,dt,rep2,rep,n0p,wall_type,wall_p_type,pmin)
int i;
int typep[NN];
int neigh[NN][NEIGHBOR];
double n[NN];
double rho[FLUID];
double x[DIM][NN],dv[DIM][NN];
double dp[NN];
double dt,rep2,rep,n0p;
int wall_type,wall_p_type;
double pmin[NN];
{
    int j,l;
    double xx,yy,rr,rr2,ddv,ddv0,ddv1;

    ddv=0.0;
    for(l=1; l<=neigh[i][0]; l++)
    {
        j=neigh[i][l];
        if(typep[j]==GHOST || typep[j]==wall_type )continue;

        xx=(x[0][i]-x[0][j]);
        yy=(x[1][i]-x[1][j]);
        rr2=xx*xx+yy*yy;
        rr=sqrt(rr2);

        ddv=dt*(dp[j]-pmin[i])/rr*weight(rr,rep)/rho2pg(rho[typep[i]],rho[typep[j]])/n0p*DIM;

        ddv0=ddv*xx/rr;
        ddv1=ddv*yy/rr;

        dv[0][i] += ddv0;
        dv[1][i] += ddv1;
    }
}

/* collision routine */
/* difference of the mass is considered */
void
collision(nump,typep,neigh,x,v,rho,r2lim,dt,rep2,col_rat,wall_type,wall_p_type)
int nump;
int typep[NN];
int neigh[NN][NEIGHBOR];
double x[DIM][NN],v[DIM][NN];
double rho[FLUID];
double r2lim,dt,rep2;
double col_rat;
int wall_type,wall_p_type;
{
    int i,j,l;
    double xx[2],rr,rr2;
    double m1,m2,mm;
    double vg[DIM],vr[DIM],vabs,vm[DIM],vrat;

    for(i=0; i<nump; i++)
    {
        if(typep[i]==GHOST)continue;

        m1=rho[typep[i]];
        for(l=1; l<=neigh[i][0]; l++)
        {
            j=neigh[i][l];
            if(j<=i)continue;
            if(typep[j]==GHOST)continue;

            xx[0]=(x[0][j]-x[0][i]);
            xx[1]=(x[1][j]-x[1][i]);
            rr2=xx[0]*xx[0]+xx[1]*xx[1];
            if(rr2<r2lim)
            {
                rr=sqrt(rr2);
                m2=rho[typep[j]];
                mm=m1+m2;
                vg[0]=(m1*v[0][i]+m2*v[0][j])/mm;
                vg[1]=(m1*v[1][i]+m2*v[1][j])/mm;
                vr[0]=m1*(v[0][i]-vg[0]);
                vr[1]=m1*(v[1][i]-vg[1]);
                vabs=(vr[0]*xx[0]+vr[1]*xx[1])/rr;

                if(vabs<0.0)continue;

                vrat=1.0+col_rat;
                vm[0]=vrat*vabs*xx[0]/rr;
                vm[1]=vrat*vabs*xx[1]/rr;
                //fprintf(tela,"particles collision: i,j,rr,v = %d,%d,%f,%f\n",i,j,rr,vabs);
                if ( typep[i]!=wall_type && typep[i]!=wall_p_type )
                {
                    v[0][i]-=vm[0]/m1;
                    v[1][i]-=vm[1]/m1;
                    x[0][i]-=dt*vm[0]/m1;
                    x[1][i]-=dt*vm[1]/m1;
                }
                if ( typep[j]!=wall_type && typep[j]!=wall_p_type )
                {
                    v[0][j]+=vm[0]/m2;
                    v[1][j]+=vm[1]/m2;
                    x[0][j]+=dt*vm[0]/m2;
                    x[1][j]+=dt*vm[1]/m2;
                }
            }
        }
    }
}

/* source term setting of the pressure Poisson equation */
void
set_source(i,typep,source,n,n0p,rho,dt)
int i;
int typep[NN];
double source[NN],n[NN];
double n0p,rho[],dt;
{
    source[i] = 1.0/dt/dt*(n[i]-n0p)/n0p;
}

/* Laplacian matrix setting for the pressure Poisson equation */
void
set_matrix(i,typep,bcon,neigh,poiss,x,n0,re,re2,lambda,rho,alpha,dt)
int i;
int typep[NN];
int bcon[NN];
int neigh[NN][NEIGHBOR];
double poiss[NN][NEIGHBOR];
double x[DIM][NN];
double n0,re,re2,lambda;
double rho[FLUID];
double alpha[FLUID],dt;
{
    FILE*matriz_=fopen("matriz_bizarra.txt","w" );
    int j,l;
    double xx,yy,rr,rr2,val;

    if(bcon[i]!=0)return;

    poiss[i][0]=0.0;
    for(l=1; l<=neigh[i][0]; l++)
    {
        j=neigh[i][l];
        if( bcon[j] == -1 )poiss[i][l]=0.0;
        else
        {

            xx=(x[0][j]-x[0][i]);
            yy=(x[1][j]-x[1][i]);
            rr2=xx*xx+yy*yy;
            rr=sqrt(rr2);

            val=2.0*DIM/lambda*weight(rr,re)/n0
                /rho2sm(rho[typep[i]],rho[typep[j]]);

            poiss[i][l] = -val;
            poiss[i][0]+=  val;
        }
    }

    poiss[i][0] += alpha[typep[i]]/dt/dt;

    int cont1,cont2;
    for(cont1 =0; cont1 < NN;  cont1++)
    {
        for(cont2 = 0; cont2 < NEIGHBOR; cont2++)
        {
            fprintf(matriz_,"%lf ",(double)poiss[cont1][cont2]);
        }
        fprintf(matriz_,"\n");
    }
    fprintf(matriz_,"\n\n\n");
}

/* incomplete Cholesky decomposition */
void
incom_decomp(nump,bcon,neigh,poiss,ic)
int nump;
int bcon[NN];
int neigh[NN][NEIGHBOR];
double poiss[NN][NEIGHBOR];
double ic[NN][NEIGHBOR];
{
    int i,j,l,mi,mj,ki,kj;
    double sum;

    for(i=0; i<nump; i++)
    {
        if(bcon[i]!=0)
            continue;
        sum=poiss[i][0];
        for(l=1; l<=neigh[i][0]; l++)
        {
            j=neigh[i][l];
            if(j>i)
                continue;
            if(bcon[j]!=0)
                continue;
            sum=sum-ic[i][l]*ic[i][l];
        }
        ic[i][0]=sqrt(sum);

        for(l=1; l<=neigh[i][0]; l++)
        {
            j=neigh[i][l];
            if(j<i)
                continue;
            if(bcon[j]!=0)
                continue;
            sum=poiss[i][l];
            for(mj=1; mj<=neigh[j][0]; mj++)
            {
                kj=neigh[j][mj];
                if(kj>=i)
                    continue;
                if(bcon[kj]!=0)
                    continue;
                for(mi=1; mi<=neigh[i][0]; mi++)
                {
                    ki=neigh[i][mi];
                    if(ki==kj)
                    {
                        sum=sum-ic[i][mi]*ic[j][mj];
                        break;
                    }
                }
            }
            ic[i][l]=sum/ic[i][0];

            for(mj=1; mj<=neigh[j][0]; mj++)
            {
                kj=neigh[j][mj];
                if(i==kj)
                {
                    ic[j][mj]=ic[i][l];
                    break;
                }
            }
        }
    }
}

/* ICCG iteration solver */
int
pcg_solver(nump,bcon,neigh,poiss,x,b,ic,eps,imin,imax,dt)
int nump;
int bcon[NN];
int neigh[NN][NEIGHBOR];
double poiss[NN][NEIGHBOR];
double x[NN],b[NN];
double ic[NN][NEIGHBOR];
double eps;
int imin,imax;
double dt;
{
    int j,k;
    double r[NN],p[NN],q[NN],s[NN];
    double aa,bb,rqo,rqn,ps;

    /*  initial setting */

    set_vector_zero(nump,s);
    set_vector_zero(nump,r);
    set_vector_zero(nump,q);
    mul_matrix_vector(nump,bcon,neigh,s,poiss,x);
    sub_vectors(nump,r,b,s,1.0);
    solver_ll(nump,bcon,neigh,q,ic,r);
    copy_vectors(nump,q,p);
    mul_vectors(nump,bcon,&rqo,r,q);

    for(k=0; k<imax; k++)
    {
        mul_matrix_vector(nump,bcon,neigh,s,poiss,p);
        mul_vectors(nump,bcon,&ps,p,s);
        aa=rqo/ps;
        add_vectors(nump,x,x,p,aa);
        sub_vectors(nump,r,r,s,aa);
        solver_ll(nump,bcon,neigh,q,ic,r);
        mul_vectors(nump,bcon,&rqn,r,q);
        bb=rqn/rqo;
        rqo=rqn;
        add_vectors(nump,p,q,p,bb);

        if ( checkconv(nump,bcon,r,eps,k,imin,imax,dt,&j) == YES) return(YES);
    }
    //fprintf(tela,"iccg iteration max limit, k,j = %d,%d\n",k,j);
    return(NO);
}

void
mul_matrix_vector(nump,bcon,list,y,mat,x)
int nump;
int bcon[NN];
int list[NN][NEIGHBOR];
double y[NN];
double mat[NN][NEIGHBOR];
double x[NN];
{
    int i,j,l;

    for(i=0; i<nump; i++)
    {
        if ( bcon[i]!=0 ) continue;
        y[i]=mat[i][0]*x[i];
        for(l=1; l<=list[i][0]; l++)
        {
            j=list[i][l];
            if(bcon[j]==-1)continue;
            y[i]=y[i]+mat[i][l]*x[j];
        }
    }
}

void
mul_vectors(nump,bcon,ans,x1,x2)
int nump;
int bcon[NN];
double *ans;
double x1[NN],x2[NN];
{
    int i;

    *ans=0.0;
    for(i=0; i<nump; i++)
    {
        if ( bcon[i]==0 ) *ans+=x1[i]*x2[i];
    }
}

void
add_vectors(nump,y,x1,x2,a)
int nump;
double y[NN],x1[NN],x2[NN];
double a;
{
    int i;

    for(i=0; i<nump; i++)y[i]=x1[i]+a*x2[i];
}

void
sub_vectors(nump,y,x1,x2,a)
int nump;
double y[NN],x1[NN],x2[NN];
double a;
{
    int i;

    for(i=0; i<nump; i++)y[i]=x1[i]-a*x2[i];
}

/* solver of a triagular matrix equation; used in ICCG solver */
void
solver_ll(nump,bcon,list,y,mat,xx)
int nump;
int bcon[NN];
int list[NN][NEIGHBOR];
double y[NN];
double mat[NN][NEIGHBOR];
double xx[NN];
{
    int i,j,l;
    double x[NN];

    copy_vectors(nump,xx,x);

    /*  forward substitution */

    for(i=0; i<nump; i++)
    {
        if(bcon[i]!=0)continue;
        for(l=1; l<=list[i][0]; l++)
        {
            j=list[i][l];
            if(j>i)continue;
            if(bcon[j]!=0)continue;
            x[i]=x[i]-mat[i][l]*y[j];
        }
        y[i]=x[i]/mat[i][0];
    }

    copy_vectors(nump,y,x);

    /*  backward substitution */

    for(i=nump-1; i>=0; i--)
    {
        if(bcon[i]!=0)continue;
        for(l=1; l<=list[i][0]; l++)
        {
            j=list[i][l];
            if(j<i)continue;
            if(bcon[j]!=0)continue;
            x[i]=x[i]-mat[i][l]*y[j];
        }
        y[i]=x[i]/mat[i][0];
    }
}

/* boundary condition setting; free surface -> Dirichlet boundary condition (P=0) */
void
set_bcon(i,typep,bcon,n,n0p,cutoff,wall_type)
int i;
int typep[NN],bcon[NN];
double n[NN];
double n0p;
double cutoff;
int wall_type;
{
    if(typep[i]==GHOST || typep[i]==wall_type)bcon[i]=-1;
    else if( n[i]/n0p < cutoff )bcon[i]=1;
    else bcon[i]=0;
}

/* checked whether Dirichlet boundaries are connected */
void
check_bcon(nump,typep,bcon,neigh,check)
int nump;
int typep[NN],bcon[NN];
int neigh[NN][NEIGHBOR];
int check[NN];
{
    int i,j,l,count;

    for(i=0; i<nump; i++)
    {
        if (bcon[i]==-1) check[i]=-1;
        else if (bcon[i]==1) check[i]=1;
        else check[i]=0;
    }

    do
    {
        count=0;
        for(i=0; i<nump; i++)
        {
            if(check[i]==1)
            {
                for(l=1; l<=neigh[i][0]; l++)
                {
                    j=neigh[i][l];
                    if(check[j]==0)check[j]=1;
                }
                check[i]=2;
                count++;
            }
        }
    }
    while (count!=0);

    for(i=0; i<nump; i++)
    {
        if(check[i]==0)fprintf(tela,"Warning no Dirichlet boundary, i=%d\n",i);
    }

}

/* convergence check in the ICCG solver */
int
checkconv(nump,bcon,r,eps,k,imin,imax,dt,j)
int nump;
int bcon[NN];
double r[NN];
double eps;
int k,imin,imax;
double dt;
int *j;
{
    int i;
    double err1,err1sum;

    *j=0;
    err1sum=0.0;
    for(i=0; i<nump; i++)
    {
        if(bcon[i]!=0)continue;
        if(r[i]>0)err1=r[i]*dt*dt;
        else err1=-r[i]*dt*dt;
        if(err1>eps)
        {
            (*j)++;
        }
        err1sum=err1sum+err1;
    }

    if ( *j==0 && k>=imin )
    {
        //fprintf(tela,"checkconv k = %d\n",k);
        return(YES);
    }
    else return(NO);
}

int
outputcheck(outputtype,time_sim,output_time,output_dt,it,interval)
char outputtype[3];
double time_sim;
double *output_time;
double output_dt;
int it,interval;
{
    if (outputtype[0]=='I')
    {
        if ( (it/interval)*interval == it ) return(YES);
        else return(NO);
    }
    else if (outputtype[0]=='T')
    {
        if ( time_sim >= *output_time )
        {
            do
            {
                *output_time = *output_time + output_dt;
            }
            while ( *output_time <= time_sim );
            return(YES);
        }
        return(NO);
    }

    return(NO);
}

void
exit_job()
{
    fclose(tela);
    printf("não foi veeeeei...\n");
    exit(1);
}

/* search minimum P for the pressure gradient calculation */
void
set_pmin(i,typep,neigh,p,pmin,wall_type)
int i;
int typep[NN];
int neigh[NN][NEIGHBOR];
double p[NN],pmin[NN];
int wall_type;
{
    int j,l;

    pmin[i]=p[i];
    for(l=1; l<=neigh[i][0]; l++)
    {
        j=neigh[i][l];
        if(typep[j]==GHOST || typep[j]==wall_type)continue;
        if(p[j]<pmin[i])pmin[i]=p[j];
    }
}

/* check of velocity change calculated in the pressure gradient term */
int
rev_limiter(nump,typep,dv,dt,disa,dtra,wall_type,wall_p_type)
int nump;
int typep[NN];
double dv[DIM][NN];
double dt,disa,dtra;
int wall_type,wall_p_type;
{
    double ddv,ratio;
    int i,j;

    j=NO;
    for(i=0; i<nump; i++)
    {
        if ( typep[i]==GHOST || typep[i]==wall_type || typep[i]==wall_p_type ) continue;

        ddv=sqrt(dv[0][i]*dv[0][i]+dv[1][i]*dv[1][i]);
        ratio=ddv*dt/disa;

        if(ratio>dtra)
        {
            //fprintf(tela,"rev_limiter,  i,typep,ratio = %d %d,%f\n",i,typep[i],ratio);
            dv[0][i]=dv[0][i]*dtra/ratio;
            dv[1][i]=dv[1][i]*dtra/ratio;
            j=YES;
        }
    }
    return(j);
}

int
dt_too_small(dt,dtmax,dtminr)
double dt,dtmax,dtminr;
{
    double val;

    if( dt < dtmax*dtminr )
    {
        fprintf(tela,"dt,dtmax,dtminr=%f,%f,%f\n",dt,dtmax,dtminr);
        return(YES);
    }
    else
    {
        val=dtmax*dtminr;
        return(NO);
    }
}

double
reset_dt(dt,ratio)
double dt,ratio;
{
    dt=dt*ratio;
    return(dt);
}

double
rho2sm(a,b)
double a,b;
{
    double c;

    c = ( a + b ) / 2.0;

    return(c);
}

double
rho2pg(a,b)
double a,b;
{
    double c;

    c = a;
    return(c);
}



