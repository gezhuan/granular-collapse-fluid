/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Stokes law

//STD
#include<iostream>
#include<math.h>
#include<random>


// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/dem/domain.h>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Vec3_t                acc;
    double                 nu;
    double                  R;
    Array<Cell *>        xmin;
    Array<Cell *>        xmax;
};

void Setup(LBM::Domain & dom, void * UD)
{
   
}

void Report(LBM::Domain & dom, void * UD)
{

}

int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    String ptype;
    bool   Render = true;
    size_t nx = 100;
    size_t ny = 50;
    size_t nz = 50;
    double nu = 0.01;
    double dx = 1.0;
    double dt = 1.0;
    double Dp = 0.1;
    double R  = 10.0;   
    double ang= 0.0;
    double Tf = 40000.0;
    double N_x  =  3;
    double N_y  =  3;
    double N_z  =  3; 
    double dtOut=  0.5;
    double height= 0.7;
    {
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> Render;       infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> nu;           infile.ignore(200,'\n');
        infile >> dx;           infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> Dp;           infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> ang;          infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
        infile >> N_x;           infile.ignore(200,'\n');
        infile >> N_y;           infile.ignore(200,'\n');
        infile >> N_z;           infile.ignore(200,'\n');
        infile >> dtOut;         infile.ignore(200,'\n');
        infile >> height;         infile.ignore(200,'\n');
    }
    
    DEM::Domain d1;

    for(size_t i=0;i<N_x;i++)
    {
        for(size_t j=0;j<N_y;j++)
        {
            for(size_t k=0;k<N_z;k++)
            {
                d1.AddSphere(-1,Vec3_t(2*R+2.1*R*i+0.1*R*(rand()/RAND_MAX),2*R+2.1*R*j+0.1*R*(rand()/RAND_MAX),2*R+2.1*R*k+0.1*R*(rand()/RAND_MAX)),R,3.0);
            }
        }
    }

    //Dom.AddSphere(-1,Vec3_t(0.5*nx*dx,0.8*ny*dx,0.5*nz*dx),R,3.0);
    d1.AddPlane(-2,Vec3_t(0.5*nx*dx,0,0.5*nz*dx),0.01,nx*dx,nz*dx,1,M_PI/2.0,&OrthoSys::e0);
    d1.AddPlane(-3,Vec3_t(0.5*nx*dx,0.5*ny*dx,nz*dx),0.01,nx*dx,ny*dx,1,0,&OrthoSys::e0);
    d1.AddPlane(-4,Vec3_t(0.5*nx*dx,0.5*ny*dx,0),0.01,nx*dx,ny*dx,1,0,&OrthoSys::e0);
    d1.AddPlane(-5,Vec3_t(nx*dx,0.5*ny*dx,0.5*nz*dx),0.01,nz*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
    d1.AddPlane(-6,Vec3_t(0,0.5*ny*dx,0.5*nz*dx),0.01,nz*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
    d1.AddPlane(-7,Vec3_t(3.0,0.5*ny*dx,0.5*nz*dx),0.01,nz*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
    d1.GetParticle(-2)->FixVeloc();
    d1.GetParticle(-3)->FixVeloc();
    d1.GetParticle(-6)->FixVeloc();
    d1.GetParticle(-4)->FixVeloc();
    d1.GetParticle(-5)->FixVeloc();
    d1.GetParticle(-7)->FixVeloc();
    Dict D1;
    for(size_t i=0;i<d1.Particles.Size();i++)
    {
        D1.Set(-i-1,"Mu Kn Kt Gn",0.3,1.0e6,1.0e6,-0.2);
        d1.Particles[i]->Ff = d1.Particles[0]->Props.m*Vec3_t(0.0,-981.0,0.0);
    }

    d1.SetProps(D1);
    d1.Alpha = 0.01;   //verlet distance

    d1.Solve(Tf,dt,dtOut,NULL,NULL,filekey.CStr(),Render,Nproc);

    size_t countdel = 0;		
		for (size_t np=0;np<d1.Particles.Size();np++)
		{
        	if (d1.Particles[np]->x(1) > height*ny*dx)
        	{
				countdel = countdel+1;
           	 	d1.Particles[np]->Tag = -8;
        	}
		}
		Array<int> delpar0;
    	Array<int> delpar1;
		Array<int> delpar2;
		Array<int> delpar3;
		Array<int> delpar4;
		Array<int> delpar5;
		
    	if (countdel > 0)
		{
    		delpar0.Push(-8);
    		d1.DelParticles(delpar0);
		}

        delpar1.Push(-2);
        d1.DelParticles(delpar1);

        // save the information from domain Dom
        d1.Save("Stage_1");
        LBM::Domain dom1(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
        dom1.Step     = 1;
        dom1.Sc       = 0.0;
        UserData dat;
        dom1.UserData = &dat;
        dat.acc      = Vec3_t(Dp,0.0,0.0);
        dat.R        = R;
        dat.nu       = nu;

        for (int i=0;i<nx;i++)
        for (int j=0;j<ny;j++)
        for (int k=0;k<nz;k++)
        {
             Vec3_t v0(0.0,0.0,0.0);
             dom1.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(1.0,v0);
        }
        dom1.Load("Stage_1");
        dom1.AddPlane(-2,Vec3_t(0.5*nx*dx,0,0.5*nz*dx),0.01,nx*dx,nz*dx,1,M_PI/2.0,&OrthoSys::e0);
        dom1.AddPlane(-3,Vec3_t(0.5*nx*dx,0.5*ny*dx,nz*dx),0.01,nx*dx,ny*dx,1,0,&OrthoSys::e0);
        dom1.AddPlane(-4,Vec3_t(0.5*nx*dx,0.5*ny*dx,0),0.01,nx*dx,ny*dx,1,0,&OrthoSys::e0);
        dom1.AddPlane(-5,Vec3_t(nx*dx,0.5*ny*dx,0.5*nz*dx),0.01,nz*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
        dom1.AddPlane(-6,Vec3_t(0,0.5*ny*dx,0.5*nz*dx),0.01,nz*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
        //dom1.AddPlane(-7,Vec3_t(3.0,0.5*ny*dx,0.5*nz*dx),0.01,nz*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
        dom1.GetParticle(-2)->FixVeloc();
        dom1.GetParticle(-3)->FixVeloc();
        dom1.GetParticle(-6)->FixVeloc();
        dom1.GetParticle(-4)->FixVeloc();
        dom1.GetParticle(-5)->FixVeloc();
        //dom1.GetParticle(-7)->FixVeloc();
        Dict B;
        for(size_t i=0;i<dom1.Particles.Size();i++)
        {
            B.Set(-i-1,"Mu Kn Kt Gn",0.3,1.0e6,1.0e6,-0.2);
            dom1.Particles[i]->Ff = dom1.Particles[0]->Props.m*Vec3_t(0.0,-981.0,0.0);
        }

        dom1.SetProps(B);
        dom1.Alpha = 0.01;   //verlet distance

        dom1.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);


    

}
MECHSYS_CATCH

