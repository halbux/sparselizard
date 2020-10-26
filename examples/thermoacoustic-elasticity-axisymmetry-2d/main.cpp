// THIS IS A RAW SPARSELIZARD CODE FILE!
// VALIDATED CODE BUT NOT CLEANED
// PLEASE REFER TO OTHER EXAMPLES FOR A SMOOTH SPARSELIZARD INTRODUCTION

// CREDITS: R. HAOUARI

// Time simulation of a thermaoacoustic wave propagation
// Cylindrical wave created by a volumic heat source ( models laser beam)
// Thermoviscouss acoustics assumed
// Mechanical coupling with a glass membrane and radiation simulated as well 


#include "sparselizardbase.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

using namespace mathop;

// Arguments are:
// r cavity, membrane radius
// depth cavity depth
// thick membrane thickness
// BLth thickness of the boundary layer

mesh createmesh(double rmb, double depth, double thick, double BLth);

void sparselizard(void)
{	

    wallclock clk;

    // Axisymmetric assumption:
    setaxisymmetry();

    //Geometry
    //---------

    // Define the geometric dimensions [m]:
    double r = 500e-6, thick = 5e-6, depth = 750e-6, BLth=10e-6;
    //axis region set to avoid weird stuff in calculation
    double caxis=1e-6;
    // The domain regions as defined in 'createmesh':
    int membrane = 1, cavity = 2, topair = 3, wall=4, laser=5, BL=6 , clamp=8,axis=9 , outair=10,outrad=11,botcoupl=12, P1=13,P2=14, topcoupl=15; ; 
    // Create the geometry and the mesh:    
    mesh mymesh = createmesh(r,depth,thick,BLth);
    // Write the mesh for display:
    mymesh.write("geometry.msh");
    // Define additional regions:
    int fluid = regionunion({cavity,outair});
    int mechacoucoupl = regionintersection({fluid,membrane});
    int slip = regionunion({botcoupl,wall});

    //normal tracking (for better BC setting... if needed)
    normal(regionunion({axis,wall,botcoupl})).write(regionunion({axis,wall,botcoupl}),"normal.vtu");
    parameter nsign;
    nsign|axis=1;
    nsign|wall=1;
    nsign|botcoupl=-1;
    expression nor=(nsign*normal(slip));
    nor.write(regionunion({axis,wall,botcoupl}),"fixed_normal.vtu");


    //Unknown fields defintion
    //--------------------------

    field u("h1xyz"),v("h1xyz"), p("h1"),T("h1"), y("y"),x ("x"),L("h1xyz");
    //Setting interpolation order :
    u.setorder(membrane, 2);
    v.setorder(cavity, 2);
    p.setorder(fluid, 1);
    T.setorder(cavity, 2);
    //Lagragian multipliers , vectorial form
    L.setorder(membrane,2);

    //gas material properties around equilibrium ( 293 K and 1 Bar ) - Linear approximation in comment
    //---------------------------------------------------------------------------------------------

    expression rho0, muB,lmb,kappa,gamma,Cp,c,alpha_p,beta,rhot;
    double R_const=8.3144598,Rs=287;
    double p0=1.013e5,T0=293;

    rho0=(p0)*0.02897/(R_const*(T0));//(p+p0)*0.02897/(R_const*(T+T0));//-(T)*0.00431073+1.20494464330308;
    lmb=-8.38278E-7+8.35717342E-8*(T0)-7.69429583E-11*pow(T0,2)+4.6437266E-14*pow(T0,3)-1.06585607E-17*pow(T0,4);//(T)*4.99221E-08+1.81322817E-05;
    muB=0.6*lmb;//(T)*2.99532E-08+1.08793690E-05;
    kappa=-0.00227583562+1.15480022E-4*(T0)-7.90252856E-8*pow(T0,2)+4.11702505E-11*pow(T0,3)-7.43864331E-15*pow(T0,4);//(T)*7.96428E-05+2.57563324E-02;
    gamma=1.4;//-(T)*1.8214E-05+1.39948988;
    Cp=1047.63657-0.372589265*(T0)+9.45304214E-4*pow(T0,2)-6.02409443E-7*pow(T0,3)+1.2858961E-10*pow(T0,4);//(T)*0.032741694+1.00541619E+03;
    c=sqrt(1.4*R_const/0.02897*(T0));//(T)*0.589990819+343.051751;
    alpha_p=sqrt(Cp*(gamma-1)/T0)/c;
    beta=gamma/(rho0*pow(c,2));
    rhot=rho0*(beta*(p+p0)-alpha_p*(T+T0));

    // Membrane properties: Young modulus, Poisson's ratio and volumic mass
    //----------------------------------------------------------------------
    //parameter E, nu, rho;
    double Es = 73.1e9,nus= 0.17;

    parameter rho;
    rho|membrane = 2203;
    rho|fluid= rho0;
    //incorporates the fluid volumic mass in the "parameter" container


    //Creation of the weak form object
    //-----------------------------------
    formulation chambre;

    //Scaling factor for numerical conditionning (if needed)
    //-----------------------------------------------------------
    double scaling = 1e0;


    // Standard isotropic elasticity
    //-------------------------------
    chambre += integral(membrane, predefinedelasticity(dof(u), tf(u),Es, nus) );
    //Inertial term:
    chambre += integral(membrane, -rho*dtdt(dof(u))*tf(u) );

    //Mass conservation
    //--------------------
    chambre += integral(cavity,(rho*div(dof(v))+rho*(beta*(dt(dof(p)))-alpha_p*(dt(dof(T)))))*tf(p) );

    //Navier-Stokes
    //--------------
    chambre += integral(cavity,-rho*dt(dof(v))*tf(v));
    chambre += integral(cavity,-grad(dof(p))*tf(v)); 
    chambre += integral(cavity,-lmb*(doubledotproduct(grad(dof(v)),grad(tf(v)))+doubledotproduct(transpose(grad(dof(v))),grad(tf(v)))));
    chambre += integral(cavity,(2*lmb/3-muB)*div(dof(v))*trace(grad(tf(v))) );
    //non linearity terms of NS
    //expression(v).reuseit();
    //chambre += integral(fluid,-rho*( grad(v)*dof(v) + grad(dof(v))*v - grad(v)*v )*tf(v));

    //thermal process
    //----------------
    chambre+=integral(cavity,kappa*grad(dof(T))*grad(tf(T)));
    chambre+=integral(cavity,(-alpha_p*T0*dt(dof(p))+rho*Cp*dt(dof(T)))*tf(T));


    //laser heat source in the cavity volume
    //-----------------------------------------
    double shift=.1e-6,pulsewidth=1e-7,waist=100e-6,absp=1,power=1;
    //spatial shape
    expression laser_shape=power*2/(getpi()*pow(waist,2))*pow(2.71828,-2*pow(x/waist,2));
    //gausian time pulse
    expression time_shape=pow(2.71828,-pow(getpi()/pulsewidth*(t()-shift),2));//
    //square time pulse
    //expression time_shape(t()-(shift-pulsewidth/2),expression(-t()+(shift+pulsewidth/2),power/pulsewidth,0),0);
    //Beert-Lambert law for absorption
    chambre+=integral(cavity,-absp*pow(2.71828,absp*y)*laser_shape*time_shape*tf(T));


    // The elastic wave propagation equation
    //------------------------------------------
    chambre += integral(outair, -grad(dof(p))*grad(tf(p)) -1/pow(c,2)*dtdt(dof(p))*tf(p));
    // A Sommerfeld condition is used on the fluid boundary to have outgoing waves:
    chambre += integral(outrad, -1/c*dt(dof(p))*tf(p));

    // Elastoacoustic coupling terms.
    //-----------------------------------
    //stress transmission from fluid to solid along the normal (carreful normal direction for the sign !!!)
    chambre += integral(botcoupl, dof(p)*normal(mechacoucoupl)*tf(u) );
    //Use of Lagrange multipliers Lx and Ly to link membrane and fluid velocity
    //Vectorial form for ease of implementation
    chambre += integral(botcoupl, -dof(L)*tf(u));
    chambre += integral(botcoupl, dof(L)*tf(v));
    chambre += integral(botcoupl, (dt(dof(u))-dof(v))*tf(L));

    //Coupling for the upper region (pure mechanical)
    chambre += integral(topcoupl, -dof(p)*normal(topcoupl)*tf(u) * scaling);
    chambre += integral(topcoupl, rho0*normal(topcoupl)*dtdt(dof(u))*tf(p)/scaling);

    //Constrained type boundary conditions
    //--------------------------------------
    //membrane clamped on the clamp line:
    u.setconstraint(clamp);
    //no slip on the wall
    v.setconstraint(wall);
    //isothermal on the interfaces
    T.setconstraint(regionunion({botcoupl,wall}));
    //symmetry (axial case : for pure 2D, not in axisym)
    //v.compx().setconstraint(axis);
    //u.compx().setconstraint(axis);
    //No Lz Lagrange multiplier
    L.compz().setconstraint(mechacoucoupl);


    //Slip BC
    //------------
    //


    //set non-zero initial values : T0 and p0
    //-----------------------------------------
    //T.setvalue(fluid,T0);
    //p.setvalue(fluid,p0);
    //init.setdata(fluid,T);
    //init.setdata(fluid,p);


    // Define the Gen alpha object to time-solve formulation 'chambre' 
    //------------------------------------------------------------------

    genalpha gena(chambre, vec(chambre), vec(chambre), 3, {false,true,true,true});
    gena.setparameter(0.75);

    gena.settolerance(1e-10);

    //variable time stepping resolution
    //--------------------------------------
    gena.setadaptivity(0.01, 10e-9, 5e-7);


    while (gettime() < 99.99e-6)
    {
        gena.next(-1);

        std::string num_file = std::to_string(gena.count());

        // Write with an order 2 interpolation and with the name of your choice (the name is also the relative path)
        u.write(membrane, "u"+num_file+".vtu",2); 
        //dtu.write(membrane, "dtu"+num_file+".vtu",2); 
        v.write(cavity, "v"+num_file+".vtu",2);		
        p.write(cavity, "pch"+num_file+".vtu",2);
        p.write(outair, "pout"+num_file+".vtu",2);
        T.write(cavity, "T"+num_file+".vtu",2); 
    }
    
    // Write .pvd files:
    grouptimesteps("u.pvd", "u", 1, gena.gettimes());
    grouptimesteps("v.pvd", "v", 1, gena.gettimes());
    grouptimesteps("pch.pvd", "pch", 1, gena.gettimes());
    grouptimesteps("pout.pvd", "pout", 1, gena.gettimes());
    grouptimesteps("T.pvd", "T", 1, gena.gettimes());

}

// THE MESH BELOW IS FULLY STRUCTURED AND IS CREATED USING THE (BASIC) SPARSELIZARD GEOMETRY CREATION TOOL.
// BOUNDARY LAYER IS ALSO INCLUDED AS WELL AS A CLOSE LINE TO THE Z-AXIS TO ACT AS THE REVOLUTION AXE

mesh createmesh(double rmb, double depth, double thick, double BLth)
{
    //closeness of the axis
    double caxis=1e-6, rarc=2e-3;
    // Give names to the physical region numbers:
    int membrane = 1, cavity = 2, topair = 3, wall=4, laser=5, BL=6 , clamp=8,axis=9 , outair=10,outrad=11,botcoupl=12, P1=13,P2=14, topcoupl=15; ; 

    // Number of mesh layers:
    int nxmb = 49, nzmb = 4, nzcavity = 73, nair = 50, nBLth=10;

    // define a vector of nodes along a vertical line

    std::vector<double>cavity_vpoints(3*(nzcavity+2*nBLth+1));

    double step1=sqrt(BLth)/nBLth, step2=(depth-2*BLth)/nzcavity;

    // nodes are spaced accordignly to have a boundary layer mesh

    for( int i=0;i<cavity_vpoints.size()/3;i++)
    {
        if(i<=nBLth)
        {
            // if linear BL :/*i*BLth/nBLth*/
            cavity_vpoints[3*i]=0;cavity_vpoints[3*i+1]=-std::pow(i*step1,2);cavity_vpoints[3*i+2]=0;	
            //std::cout<<i<<" "<<-std::pow(i*step1,2)<<std::endl;
        }
        else
        {
            if(i> nBLth && i<=nBLth+nzcavity)
            {
                cavity_vpoints[3*i]=0;cavity_vpoints[3*i+1]=-BLth-(i-nBLth)*step2;cavity_vpoints[3*i+2]=0;	
            }
            else
            {
                // if linear BL :-(i-nzcavity-2*nBLth)*BLth/nBLth
                cavity_vpoints[3*i]=0;cavity_vpoints[3*i+1]=-depth+std::pow((i-nzcavity-2*nBLth)*step1,2);cavity_vpoints[3*i+2]=0;		
            }
        }

        }
    cavity_vpoints[3*(nzcavity+2*nBLth)]=0;cavity_vpoints[3*(nzcavity+2*nBLth)+1]=-depth;cavity_vpoints[3*(nzcavity+2*nBLth)+2]=0;	

    // define a vector of nodes along a horizontal line

    std::vector<double>cavity_hpoints(3*(2+nxmb+nBLth+1));

    step1=sqrt(BLth)/nBLth;
    step2=(rmb-BLth-caxis)/(nxmb);

    //special node for the axis ; should be smaller than one node

    cavity_hpoints[0]=0;cavity_hpoints[1]=0;cavity_hpoints[2]=0;
    cavity_hpoints[3]=caxis;cavity_hpoints[4]=0;cavity_hpoints[5]=0;

    for( int i=2;i<cavity_hpoints.size()/3;i++)
    {
        if(i<=nxmb)
        {
            cavity_hpoints[3*i]=caxis+(i-1)*step2;cavity_hpoints[3*i+1]=0;cavity_hpoints[3*i+2]=0;			
        }
        else
        {
            // if linear BL : (i-2-nxmb-nBLth)*BLth/nBLth
            cavity_hpoints[3*i]=rmb-std::pow((i-2-nxmb-nBLth)*step1,2);cavity_hpoints[3*i+1]=0;cavity_hpoints[3*i+2]=0;	
        }

    }

    // define lines with nodes predefined in previous vectors

    shape ligne1("line",axis,cavity_vpoints);
    shape ligne2("line",botcoupl,cavity_hpoints);
    shape ligne3("line",-1,cavity_vpoints);
    shape ligne4("line",-1,cavity_hpoints);

    ligne3.shift(rmb,0,0);
    ligne4.shift(0,-depth,0);
    //printvector(ligne1.getcoords());
    //printvector(ligne2.getcoords());

    // derive rectangle with boundary layer

    shape nK("quadrangle",cavity,{ligne1,ligne4,ligne3,ligne2});
    //shape nK("quadrangle",cavity,{0,0,0,0,-depth,0,rmb,-depth,0,rmb,0,0},{nzcavity,nxmb,nzcavity,nxmb});

    shape ligne5("line",topcoupl,cavity_hpoints);
    ligne5.shift(0,thick,0);
    shape ligne6("line",clamp,{rmb,0,0,rmb,thick,0},nzmb);
    shape ligne7("line",-1,{rmb,0,0,rmb,thick,0},nzmb);
    ligne7.shift(-rmb,0,0);

    // axis line defined for formulation

    shape mb("quadrangle",membrane,{ligne7,ligne2,ligne6,ligne5});
    //shape laxe=ligne1.duplicate();
    //laxe.shift(caxis,0,0);
    //laxe.setphysicalregion(axis);

    // walls definition

    shape walline("union",wall,{ligne4,ligne3});

    //std::cout<<cavity_hpoints.size()/3<<"   "<< 3+nxmb+nBLth<<"    "<<sin(3.14/4)<<endl;

    double angle=3.1415/20;

    shape quart("arc",outrad,{rarc*cos(angle),rarc*sin(angle)+thick,0,0,thick+rarc,0,0,thick,0},3+nxmb+nBLth);
    shape ligne8("line",-1,{0,thick,0,0,rarc,0},nair);
    shape ligne9("line",outrad,{rmb,thick,0,rarc*cos(angle),thick+rarc*sin(angle),0},nair);

    shape out("quadrangle",outair,{ligne8,ligne5,ligne9,quart});

    shape Po1 = ligne2.getsons()[0];
    Po1.setphysicalregion(P1);

    shape Po2 = ligne2.getsons()[1];
    Po2.setphysicalregion(P2);

    mesh mymesh({nK,mb,ligne1,ligne2,ligne5,ligne6,walline,quart,ligne9,out,Po1,Po2});
    //mymesh.write("test.msh");


    return mymesh;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

