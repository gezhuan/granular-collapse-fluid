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

// MechSys
#include <mechsys/lbm/Domain.h>

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
    /*double delta;
    delta = dom.Particles[i]->R - dom.Particles[i]->X(1);
    if (delta > 0.0)
    {
        dom.Particles[i]->Ff(1) += dat.Kn*delta-dat.Gn*dom.Particles[i]->V(1);
        dat.Fb[i] += dom.Particles[i]->V(0)*dom.dt;
        if (fabs(dat.Fb[i])>=dom.Particles[i]->Mu*delta) dat.Fb[i] *= dom.Particles[i]->Mu*delta/fabs(dat.Fb[i]);
        dom.Particles[i]->Ff(0) += -dat.Kn*dat.Fb[i];
    }*/
    /*UserData & dat = (*static_cast<UserData *>(UD));
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c   = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.acc;
    }
    for (size_t i=0;i<dat.xmin.Size();i++)
    {
        Cell * c = dat.xmin[i];
        if(c->IsSolid) continue;
        c->F[1] = 1.0/3.0 *(-2*c->F[0] - 4*c->F[10] - 4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[7] = 1.0/24.0*(-2*c->F[0] - 4*c->F[10] - 4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
        c->F[9] = 1.0/24.0*(-2*c->F[0] + 20*c->F[10] - 4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[11]= 1.0/24.0*(-2*c->F[0] - 4*c->F[10] + 20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[13]= 1.0/24.0*(-2*c->F[0] - 4*c->F[10] - 4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.xmax.Size();i++)
    {
        Cell * c = dat.xmax[i];
        if(c->IsSolid) continue;
        c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
        c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
        c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }*/
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
    }
    
    

    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    Dom.Step     = 1;
    Dom.Sc       = 0.0;
    UserData dat;
    Dom.UserData = &dat;
    dat.acc      = Vec3_t(Dp,0.0,0.0);
    dat.R        = R;
    dat.nu       = nu;

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(1.0,v0);
    }  

    for(size_t i=0;i<N_x;i++)
    {
        for(size_t j=0;j<N_y;j++)
        {
            for(size_t k=0;k<N_z;k++)
            {
                Dom.AddSphere(-i-1,Vec3_t(2+2*R*i+0.5*R*(rand()/RAND_MAX),2+2*R*j++0.5*R*(rand()/RAND_MAX),2+2*R*k++0.5*R*(rand()/RAND_MAX)),R,3.0);
            }
        }
    }

    //Dom.AddSphere(-1,Vec3_t(0.5*nx*dx,0.8*ny*dx,0.5*nz*dx),R,3.0);
    Dom.AddPlane(-N_x*N_y*N_z-1,Vec3_t(0.5*nx*dx,0,0.5*nz*dx),1,nx*dx,nz*dx,1,M_PI/2.0,&OrthoSys::e0);
    Dom.AddPlane(-N_x*N_y*N_z-2,Vec3_t(0.5*nx*dx,0,nz*dx),1,nx*dx,ny*dx,1,0,&OrthoSys::e0);
    Dom.AddPlane(-N_x*N_y*N_z-3,Vec3_t(0.5*nx*dx,0,0),1,nx*dx,ny*dx,1,0,&OrthoSys::e0);
    Dom.AddPlane(-N_x*N_y*N_z-4,Vec3_t(0,0,0.5*nx*dx),1,nx*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
    Dom.AddPlane(-N_x*N_y*N_z-5,Vec3_t(nx*dx,0,0.5*nx*dx),1,nx*dx,ny*dx,1,M_PI/2.0,&OrthoSys::e1);
    Dom.GetParticle(-N_x*N_y*N_z-2)->FixVeloc();
    Dom.GetParticle(-N_x*N_y*N_z-3)->FixVeloc();
    Dom.GetParticle(-N_x*N_y*N_z-1)->FixVeloc();
    Dom.GetParticle(-N_x*N_y*N_z-4)->FixVeloc();
    Dom.GetParticle(-N_x*N_y*N_z-5)->FixVeloc();
    Dict D;
    for(size_t i=0;i<Dom.Particles.Size();i++)
    {
        D.Set(-i-1,"Mu Kn Kt Gn",0.3,1.0e6,1.0e6,-0.2);
        Dom.Particles[i]->Ff = Dom.Particles[0]->Props.m*Vec3_t(0.0,-981.0,0.0);
    }
    //D.Set(-N_x*N_y*N_z-1,"Mu Kn Kt Gn",0.3,1.0e6,1.0e6,-0.2);
    //D.Set(-N_x*N_y*N_z-2,"Mu Kn Kt Gn",0.3,1.0e6,1.0e6,-0.2);
    //D.Set(-N_x*N_y*N_z-3,"Mu Kn Kt Gn",0.3,1.0e6,1.0e6,-0.2);
    Dom.SetProps(D);
    Dom.Alpha = 0.01;   //verlet distance

    /*for (size_t i=0;i<nx;i++)
    {
        for(size_t j=0;j<nz;j++)
        {
            Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
            Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
        }      
    }*/
    //Dom.Particles[0]->FixVeloc();
    //Dom.Particles[0]->w = Vec3_t(0.0,0.0,w);
    //Dom.Particles[0]->Ff = Dom.Particles[0]->Props.m*Vec3_t(0.0,-981.0,0.0);
    //Solving
    Dom.Solve(Tf,0.01*Tf,Setup,Report,filekey.CStr(),Render,Nproc);
}
MECHSYS_CATCH

