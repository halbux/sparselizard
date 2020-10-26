//2D fluid convection simulation
//Heat source : 2D projection of a hand
//fluid is considered compressible
//hence density is not constant and div !=0
// 
// WORK IN PROGRESS, NOT VALIDATED YET! Help is appreciated!
//
// CREDITS. R. Haouari, A. Halbach


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
	//Mesh import
	mesh mymesh("natconv.msh");
	
	// Init state : 1atm and 20°C  Hand temperature :37°C
	double p0=1.01325e5,T0=293.15,Tc=(37+273.15-T0);

	// The domain regions as defined in 'hand_geom.geo':
    int fluid = 1, disk = 2, inlet = 3, outlet = 4, sides = 5;
	
	int all = regionunion({fluid, disk});
	int diskskin = regionintersection({fluid, disk});
	
	int zflux = regionunion({sides, outlet});
	int noslip = regionunion({sides, diskskin});
	int fluidskin = regionunion({sides, diskskin, inlet, outlet});

    // Field v is the flow velocity. It uses nodal shape functions "h1" with two components in 2D.
    // Field p and T are the relative pressure andtemperature respectively. 
	// x et y speak by themselves 
    field v("h1xy"), p("h1"), T("h1"), x("x"),y("y");
	
	// Use order 2 for better precision
	//LBB condition overidden with the use od P+P stabilization
    p.setorder(all, 2); 
	v.setorder(all, 2); 
    T.setorder(all, 2);
	
		
	//normal tracking and fixing
	normal(fluidskin).write(fluidskin, "normal.pos", 2);
		
	//fluid material properties around equilibrium ( 293 K and 1 Bar )
	//----------------------------------------------------------------

	expression rho0, muB,lmb,kappa,Cp,rhot,dtrho,gradrho;
	double R_const=8.3144598,Rs=287;
		
	lmb=-8.38278E-7+8.35717342E-8*(T0)-7.69429583E-11*pow(T0,2)+4.6437266E-14*pow(T0,3)-1.06585607E-17*pow(T0,4);//(T)*4.99221E-08+1.81322817E-05;
	muB=0.6*lmb;
	kappa=-0.00227583562+1.15480022E-4*(T0)-7.90252856E-8*pow(T0,2)+4.11702505E-11*pow(T0,3)-7.43864331E-15*pow(T0,4);//(T)*7.96428E-05+2.57563324E-02;
	Cp=1047.63657-0.372589265*(T0)+9.45304214E-4*pow(T0,2)-6.02409443E-7*pow(T0,3)+1.2858961E-10*pow(T0,4);//(T)*0.032741694+1.00541619E+03;
	rho0=(p0)*0.02897/(R_const*(T0));
	rhot=(p0+p)*0.02897/(R_const*(T+T0));
	dtrho=rho0*(dt(p)/(p+p0)-dt(T)/(T+T0));
	gradrho=rho0*(grad(p)/(p+p0)-grad(T)/(T+T0));
	
	//Dirichlet BC
	//-----------------
	// no-slip (0 velocity) condition 
    v.setconstraint(noslip);	
	//in and outflow
	v.setconstraint(inlet,1e-3*array2x1(0,1));
	p.setconstraint(outlet);
	
	//warm disk
	T.setconstraint(disk, Tc);
	//fresh air in
	T.setconstraint(inlet);
	
	

	// Hp-adaptivity:
	expression criterion = zienkiewiczzhu(grad(compx(v))) + zienkiewiczzhu(grad(compy(v)));
	
	double errtarget = 1.0;
	expression errorindicator = ifpositive(criterion - errtarget, 1, 0);
	
	T.setvalue(disk, Tc);
	
	expression crit = ifpositive( abs(T) - 1.0, 1, 0 );
	
	mymesh.setadaptivity(crit, 0, 3);
	
	// Make mesh refined around disk for first iteration:
	for (int i = 0; i < 3; i++)
	    adapt();
	
	crit = errorindicator;
	
	mymesh.setadaptivity(crit, 0, 3);
	
	v.setorder(crit, 2, 3); // FIXME: Even with p-stab I cannot use order 1 for v. Why?
	p.setorder(crit, 1, 2);
	T.setorder(crit, 1, 2);
	


    // Define the weak formulation for the convection:
    formulation convec;
	
    //flow
	//-----
    // Define the weak formulation for time-dependent compressible laminar flow:
    convec += integral(fluid, predefinednavierstokes(dof(v), tf(v),v, dof(p), tf(p), lmb, rhot, dtrho, gradrho, true,false,true)) ;
	//gravitational force : responsible for the bouyancy - convection
	convec += integral(fluid, (rhot)*9.81*array2x1(0,-1)*tf(v));
	//P+P stab
	convec+=integral(fluid,1e-10*grad(dof(p))*grad(tf(p)));
		
	
	//thermal process
	//----------------
	//conduction
	convec+=integral(fluid,kappa*grad(dof(T))*grad(tf(T)));
	//time variation
	convec+=integral(fluid,rho0*Cp*dt(dof(T))*tf(T));
	//advection part 1
	convec+=integral(fluid,rho0*Cp*v*grad(dof(T))*tf(T));
	//advection part 2 - compressible - can be removed if fluid incompressible
	convec+=integral(fluid,rho0*Cp*div(v)*dof(T)*tf(T));
	//zero flux BC
	convec+=integral(fluid,-(-kappa*grad(dof(T)))* normal(fluidskin) *tf(T,zflux));  // WHAT? normal --> elements devraient etre une boundary! pas fluid FIXME
	    
    // Initial state: Tc on disk, 0 on fluid and 0 dt(T):
	T.setvalue(fluid);
	vec dtinit(convec);
	
	// implicit Euler for time resolution 
    // initiation with the dtinit vector
    impliciteuler eul(convec, dtinit);
	//set solver parameters
	eul.settolerance(1e-3);
	eul.setrelaxationfactor(.55);
	
	// Time step and number of time steps:
	double ts = 0.01;
	int nts = 1000;
			
	settime(0);
	for (int i = 0; i < nts; i++)
	{
		// Call the nonlinear solver and save the vector of solution
		eul.next(ts, 20);

		v.write(fluid, "v"+std::to_string(1000 + i)+".vtu", 3);		
		p.write(fluid, "p"+std::to_string(1000 + i)+".vtu", 3);
		T.write(fluid, "T"+std::to_string(1000 + i)+".vtu", 3);
		
		criterion.write(fluid, "criterion"+std::to_string(1000 + i)+".vtu", 2);
		fieldorder(v.compx()).write(fluid, "fieldorderv"+std::to_string(1000 + i)+".vtu", 2);
	    
	    double maxcrit = criterion.max(fluid, 1)[0];
	    
		std::cout << std::endl << "@" << i << " -----------------------------------------------> " << maxcrit << " (#dofs is " << convec.countdofs() << ")" << std::endl;
		
	    adapt(1);
	}
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

