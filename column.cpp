/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
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

#include <math.h>
#include <random>

#include <gsl/gsl_linalg.h>
// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

struct UserData
{
    double visc    = 1.0e-4;  //viscosity of fluid, what is roughly 1e-7
    double hcutoff = 0.25;  // h cut off
	double hrough  = 0.01;  // surface roughness
};

void Setup (DEM::Domain & dom, void * UD) //function that is going to be called every time step
{
    UserData & dat = (*static_cast<UserData *>(UD));
    
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0; i<dom.Interactons.Size(); i++)
    {
        DEM::Particle * P1 = dom.Interactons[i]->P1;
        DEM::Particle * P2 = dom.Interactons[i]->P2;
        double R1 = P1->Props.R;
        double R2 = P2->Props.R;
		double Reff = R1*R2/(R1+R2);

        double dist = norm(P1->x - P2->x) - R1 - R2;
		
        if (dist > dat.hcutoff) continue;
		if (dist < dat.hrough) continue;
        Vec3_t normal = P2->x - P1->x;
        normal /= norm(normal);

        Vec3_t t1,t2,x1,x2;
        Rotation(P1->w,P1->Q,t1);
        Rotation(P2->w,P2->Q,t2);
        x1 = R1*normal;
        x2 = -R2*normal;
        Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
        Vec3_t vn = dot(normal,vrel)*normal;
        Vec3_t vt = vrel - vn;

        Vec3_t Fn = -6.0*3.1415926*Reff*Reff*dat.visc*vn/dist;
        Vec3_t Ft = -6.0*3.1415926*Reff*dat.visc*log(8.0/15.0*(Reff/dist)+0.9588)*vt; // /norm(vt);
        t1 = cross(x1,Ft);
        t2 = -cross(x2,Ft);
        
        Vec3_t t1p,t2p;
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (t1,q,t1p);
        Conjugate (P2->Q,q);
        Rotation  (t2,q,t2p);


        omp_set_lock  (&dom.Interactons[i]->P1->lck);
        dom.Interactons[i]->P1->F += Fn + Ft;
        dom.Interactons[i]->P1->T += t1p;
        omp_unset_lock(&dom.Interactons[i]->P1->lck);
        omp_set_lock  (&dom.Interactons[i]->P2->lck);
        dom.Interactons[i]->P2->F += -Fn - Ft;
        dom.Interactons[i]->P2->T += t2p;
        omp_unset_lock(&dom.Interactons[i]->P2->lck);
    }

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0; i<dom.Particles.Size(); i++)
    {
		
        dom.Particles[i]->F += -0.0*6.0*3.1415926*(dat.hcutoff/2.0)*dat.visc*dom.Particles[i]->v;
    }
}


int main(int argc, char **argv) try
{   
	if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
	
	// Setting number of CPUs
	size_t Nproc = 1;
	if (argc>=3) Nproc = atoi(argv[2]);
	
	String filekey  (argv[1]);
	String filename (filekey+".inp");
	if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
	ifstream infile(filename.CStr());
    
    // domain
    UserData dat;
    DEM::Domain d;
    dat.visc = 1.0e-5;
	
	String CrossSection;// Shape of the cross-section of the column
	String ptype;       // Particle type 
	String test;       // Test type 
	bool   Cohesion;    // Decide if coheison is going to be simulated
	double fraction;    // Fraction of particles to be generated
	double Kn;          // Normal stiffness
	double Kt;          // Tangential stiffness
	double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Muw;         // Frictional coefficient of the bottom wall
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
	double etaf;        // dynamics viscosity of interstitial fluid
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
	double dtOut1;      // Time step for output for the dropping stage
    double dtOut;       // Time step for output for the collapsing stage
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
	double Lx_p;         // 
    size_t scalingx;    // scalingx
    size_t scalingy;    // scalingy
    size_t scalingz;    // scalingz
    size_t plane_x;     // the scaling of the size of the plane in x direction, how many particles per unit length
    size_t plane_y;     // the scaling of the size of the plane in y direction, how many particles per unit length
    double rho;         // rho
	double Tf1;         // Final time for the dropping stage test
    double Tf;          // Final time for the collapsing test
	double height;
    {
		infile >> CrossSection;     infile.ignore(200,'\n');
		infile >> ptype;     infile.ignore(200,'\n');
		infile >> test;     infile.ignore(200,'\n');
		infile >> Cohesion;     infile.ignore(200,'\n');
		infile >> fraction;     infile.ignore(200,'\n');
		infile >> Kn;           infile.ignore(200,'\n');
		infile >> Kt;           infile.ignore(200,'\n');
		infile >> Gn;           infile.ignore(200,'\n');
		infile >> Gt;           infile.ignore(200,'\n');
		infile >> Mu;           infile.ignore(200,'\n');
		infile >> Muw;          infile.ignore(200,'\n');
		infile >> Bn;           infile.ignore(200,'\n');
		infile >> Bt;           infile.ignore(200,'\n');
		infile >> Bm;           infile.ignore(200,'\n');
		infile >> Eps;          infile.ignore(200,'\n');
		infile >> etaf;         infile.ignore(200,'\n');
		infile >> R;            infile.ignore(200,'\n');
		infile >> seed;         infile.ignore(200,'\n');
		infile >> dt;           infile.ignore(200,'\n');
		infile >> dtOut1;       infile.ignore(200,'\n');
		infile >> dtOut;        infile.ignore(200,'\n');
		infile >> Lx;           infile.ignore(200,'\n');
		infile >> Ly;           infile.ignore(200,'\n');
		infile >> Lz;           infile.ignore(200,'\n');
		infile >> Lx_p;         infile.ignore(200,'\n');
		infile >> scalingx;     infile.ignore(200,'\n');
		infile >> scalingy;     infile.ignore(200,'\n');
		infile >> scalingz;     infile.ignore(200,'\n');
		infile >> plane_x;      infile.ignore(200,'\n');
		infile >> plane_y;      infile.ignore(200,'\n');
		infile >> rho;          infile.ignore(200,'\n');
		infile >> Tf1;          infile.ignore(200,'\n');
		infile >> Tf;           infile.ignore(200,'\n');
		infile >> height;       infile.ignore(200,'\n');
    }

	// dat.visc = 0.0;   //set properties of fluid
    //Some key parameters
    //size_t Nx = size_t(Lx_p*scalingx);  //Division of the rectangular box
    //size_t Ny = size_t(Ly*scalingy);
    //size_t Nz = size_t(Lz*scalingz);
	double Nx = Lx_p*scalingx-2;
	double Ny = Ly*scalingy-2;
	double Nz = Lz*scalingz*1.8;
    Kn = Kn/(scalingx*scalingy); //Stiffness constant for particles
    Kt = Kt/(scalingx*scalingy);

    //Add the granular column
    
	Vec3_t axis0(OrthoSys::e0); // rotation of face
	Vec3_t axis1(OrthoSys::e1); // rotation of face
	double Cf = 15.0;
	double plate_t = 0.02;
	
    
	// estimate the number of particles we need
	size_t num_of_particles = (scalingx*Lx_p+1)*(scalingy*Ly+1)*(scalingz*Lz*0.5);
	double Sphere_size = 1.0/scalingx;
	dat.hcutoff = Sphere_size;
	/*double delta_march = Sphere_size*1.2;
	// position of the first sphere
	Vec3_t X_Sphere (delta_march, delta_march, delta_march);
	d.AddSphere(-1, X_Sphere + Vec3_t(0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX), Sphere_size/2.0, rho);
	// go over all the particles to 
	for (size_t np=1;np<num_of_particles;np++)
	{
		if (X_Sphere(1)<Ly-delta_march)
		{
			X_Sphere(1) = X_Sphere(1) + delta_march;
			d.AddSphere(-1, X_Sphere + Vec3_t(0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX), Sphere_size/2.0, rho);
		}
		else
		{
			X_Sphere(0) = X_Sphere(0) + delta_march;
			//X_Sphere(1) = -0.5*Ly+delta_march;
			if (X_Sphere(0) < Lx_p-delta_march)
			{
				d.AddSphere(-1, X_Sphere + Vec3_t(0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX), Sphere_size/2.0, rho);
			}
			else
			{
				//X_Sphere(0) = -0.5*Lx+delta_march;
				X_Sphere(2) = X_Sphere(2)+delta_march;
				d.AddSphere(-1, X_Sphere + Vec3_t(0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX, 0.1*Sphere_size*rand()/RAND_MAX), Sphere_size/2.0, rho);
			}
		}
				
	}*/
	for(size_t i=0;i<Nx;i++)
    {
        for(size_t j=0;j<Ny;j++)
        {
            for(size_t k=0;k<Nz;k++)
            {
                d.AddSphere(-1,Vec3_t(Sphere_size+1.1*Sphere_size*i+0.1*Sphere_size*(rand()/RAND_MAX),Sphere_size+1.1*Sphere_size*j+0.1*Sphere_size*(rand()/RAND_MAX),Sphere_size+1.1*Sphere_size*k+0.1*Sphere_size*(rand()/RAND_MAX)),Sphere_size/2.0,rho);
            }
        }
    }
		
	d.AddPlane (-11, Vec3_t(Lx,Ly/2.0,Lz/2),  0.01, Lz, Ly, 1.0, M_PI/2.0, &axis1);
	d.AddPlane (-12, Vec3_t(0,Ly/2.0,Lz/2.0), 0.01, Lz, Ly, 1.0, 3.0*M_PI/2.0, &axis1);
	d.AddPlane (-13, Vec3_t(Lx/2.0,0.0,Lz/2.0),  0.01, Lx, Lz, 1.0, 3.0*M_PI/2.0, &axis0);
	d.AddPlane (-14, Vec3_t(Lx/2.0,Ly,Lz/2.0), 0.01, Lx, Lz, 1.0, M_PI/2.0, &axis0);
	d.AddPlane (-15, Vec3_t(Lx/2.0,Ly/2.0,0), 0.01, Lx, Ly, 1.0);
	d.AddPlane (-16, Vec3_t(Lx_p,Ly/2.0,Lz/2),  0.01, Lz, Ly, 1.0, M_PI/2.0, &axis1);
	d.GetParticle(-11)->FixVeloc();
	d.GetParticle(-12)->FixVeloc();
	d.GetParticle(-13)->FixVeloc();
	d.GetParticle(-14)->FixVeloc();
	d.GetParticle(-15)->FixVeloc();
	d.GetParticle(-16)->FixVeloc();
		
	for (size_t np=0;np<d.Particles.Size();np++)
    {
		// d.Particles[np]->Tag = -1;
		d.Particles[np]->Ff = d.Particles[np]->Props.m*Vec3_t(0.0,0.0,-981.0);
		d.Particles[np]->Props.Kn = Kn; // normal stiffness
		d.Particles[np]->Props.Kt = Kt; // trangential stiffness
		d.Particles[np]->Props.Gn = Gn; // restitution coefficient
		d.Particles[np]->Props.Mu = Mu; // frictional coefficient
	}
        
	// solve to get random packed cubes
	dt = 0.5*d.CriticalDt(); //Calculating time step
	d.Alpha = 0.05; //Verlet distance
	//d.WriteXDMF("test");
	d.Solve(/*tf*/Tf1, dt, /*dtOut*/dtOut1, NULL, NULL, "drop_spheres", 2, Nproc);

		
	size_t countdel = 0;
	Array<int> delpar0;
    Array<int> delpar1;
	Array<int> delpar2;
	Array<int> delpar3;
	Array<int> delpar4;
	Array<int> delpar5;
	Array<int> delpar6;

	delpar1.Push(-11);
	d.DelParticles(delpar1);
	d.WriteXDMF("test1");

	delpar2.Push(-12);
	d.DelParticles(delpar2);
	d.WriteXDMF("test2");

	delpar3.Push(-13);
	d.DelParticles(delpar3);
	d.WriteXDMF("test3");

	delpar4.Push(-14);
	d.DelParticles(delpar4);
	d.WriteXDMF("test4");

	delpar5.Push(-15);
	d.DelParticles(delpar5);
	d.WriteXDMF("test5");

	delpar6.Push(-16);
	d.DelParticles(delpar6);
	d.WriteXDMF("test6");

	for (size_t np=0;np<d.Particles.Size();np++)
	{
        if (d.Particles[np]->x(0) > Lx || d.Particles[np]->x(0) < 0 || d.Particles[np]->x(1) > Ly || d.Particles[np]->x(1) < 0 || d.Particles[np]->x(2) > Lz*height || d.Particles[np]->x(2) < 0)
        {
			countdel = countdel+1;
           	 d.Particles[np]->Tag = 10;
        }
	}

		
    if (countdel > 0)
	{
    	delpar0.Push(10);
    	d.DelParticles(delpar0);
	}


	// save the information from domain d
	d.Save("Stage_1");


    // solve the problem, but there's something different for the cubes   
    if (ptype=="sphere" || ptype=="Sphere")
    {
		
        DEM::Domain dom(&dat);
		dat.visc = etaf; // set properties of fluid
        dom.Load("Stage_1");
        //Determinaning the bounding box
        Vec3_t Xmin,Xmax;
        dom.BoundingBox(Xmin,Xmax);

        //Adding plate at the base of the column
        //dom.AddPlane(-2,Vec3_t(0.0,0.0,Xmin(2)-plate_t),plate_t,plane_x*Lz,plane_y*Lz,rho);
		dom.AddPlane (-11, Vec3_t(Lx,Ly/2.0,Lz/2),  0.01, Lz, Ly, 1.0, M_PI/2.0, &axis1);
	    dom.AddPlane (-12, Vec3_t(0,Ly/2.0,Lz/2.0), 0.01, Lz, Ly, 1.0, 3.0*M_PI/2.0, &axis1);
	    dom.AddPlane (-13, Vec3_t(Lx/2.0,0.0,Lz/2.0),  0.01, Lx, Lz, 1.0, 3.0*M_PI/2.0, &axis0);
	    dom.AddPlane (-14, Vec3_t(Lx/2.0,Ly,Lz/2.0), 0.01, Lx, Lz, 1.0, M_PI/2.0, &axis0);
	    dom.AddPlane (-15, Vec3_t(Lx/2.0,Ly/2.0,0), 0.01, Lx, Ly, 1.0);
	    
	    dom.GetParticle(-11)->FixVeloc();
	    dom.GetParticle(-12)->FixVeloc();
	    dom.GetParticle(-13)->FixVeloc();
	    dom.GetParticle(-14)->FixVeloc();
	    dom.GetParticle(-15)->FixVeloc();
	    

        //Fixing the Plane so it does not move (plane tag is -2)
        //dom.GetParticle(-2)->FixVeloc();
        // set properties for the particles
        for (size_t np=0;np<dom.Particles.Size();np++)
    	{
        	// d.Particles[np]->Tag = -1;
			dom.Particles[np]->Ff = dom.Particles[np]->Props.m*Vec3_t(0.0,0.0,-981.0);
        	dom.Particles[np]->Props.Kn = Kn; // normal stiffness
         	dom.Particles[np]->Props.Kt = Kt; // trangential stiffness
        	dom.Particles[np]->Props.Gn = Gn; // restitution coefficient
         	dom.Particles[np]->Props.Mu = Mu; // frictional coefficient
    	}
	    // set the frictional coefficient for the bottom wall
	
	    Dict B;
	
        B.Set(-2,"Mu",Muw);
	    dom.SetProps(B);
	
	    // Change the shape of cross-section
	    if (CrossSection=="circle" || CrossSection=="Circle")
	    {
		    for (size_t np=0;np<dom.Particles.Size();np++)
    	    {
        	    if (dom.Particles[np]->x(0)*dom.Particles[np]->x(0)+dom.Particles[np]->x(1)*dom.Particles[np]->x(1)>=0.25*Lx*Ly)
        	    {
           	 	    dom.Particles[np]->Tag = 10;
        	    }
    	    }
    	    Array<int> delpar;
    	    delpar.Push(10);
    	    dom.DelParticles(delpar);
	    }
	    else if (CrossSection=="right_triangle")
	    {
		    for (size_t np=0;np<dom.Particles.Size();np++)
    	    {
        	    if (dom.Particles[np]->x(1) > Ly/Lx* dom.Particles[np]->x(0))
        	{
           	 	dom.Particles[np]->Tag = 10;
        	}
    	}
    	Array<int> delpar;
    	delpar.Push(10);
    	dom.DelParticles(delpar);
	    }
	    else if (CrossSection=="isoscele_triangle")
	    {
		    for (size_t np=0;np<dom.Particles.Size();np++)
    	    {
        	    if ((dom.Particles[np]->x(1) > 2*Ly/Lx* dom.Particles[np]->x(0) + Ly/2) || (dom.Particles[np]->x(1) > -2*Ly/Lx* dom.Particles[np]->x(0) + Ly/2))
        	    {
           	    	dom.Particles[np]->Tag = 10;
        	    }
    	    }
    	    Array<int> delpar;
    	    delpar.Push(10);
    	    dom.DelParticles(delpar);
	    }
	    else if (CrossSection=="square" || CrossSection=="Square")
	    {
		    std::cout << "The cross-section is a square" << std::endl;
	    }
	    else throw new Fatal("Packing for particle type not implemented yet");

        // solve
        dt = 0.5*dom.CriticalDt(); //Calculating time step
        dom.Alpha = R; //Verlet distance
	    //d.WriteXDMF("test");
        dom.Solve(/*tf*/Tf, dt, /*dtOut*/dtOut, &Setup, NULL, "column_sphere", 2, Nproc);
    }
    else throw new Fatal("Packing for particle type not implemented yet");

        
    
}
MECHSYS_CATCH
