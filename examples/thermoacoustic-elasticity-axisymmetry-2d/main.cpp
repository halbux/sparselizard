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

    //Save folder
    std::string savefolder="./data/";


    wallclock clk;
    clock_t start=clock();

    // Axisymmetric assumption:
    setaxisymmetry();

    //Geometry
    //---------

    // Define the geometric dimensions [m]:
    double r = 500e-6, thick = 5e-6, depth = 750e-6, BLth=10e-6;
    //axis region set to avoid weird stuff in calculation
    double caxis=1e-6;
    // The domain regions as defined in 'createmesh':
    int membrane = 1, KviT = 2, topair = 3, wall=4, laser=5, BL=6 , clamp=8,axis=9 , outair=10,outrad=11,botcoupl=12, P1=13,P2=14, topcoupl=15; ; 
    // Create the geometry and the mesh:    
    mesh mymesh = createmesh(r,depth,thick,BLth);
    // Write the mesh for display:
    mymesh.write(savefolder+"geometry.msh");
    // Define additional regions:
    int fluid = regionunion({KviT,outair});
    int mechacoucoupl = regionintersection({fluid,membrane});
    int slip = regionunion({botcoupl,wall});

    //normal tracking (for better BC setting... if needed)
    normal(regionunion({axis,wall,botcoupl})).write(regionunion({axis,wall,botcoupl}),savefolder+"normal.vtu");
    parameter nsign;
    nsign|axis=1;
    nsign|wall=1;
    nsign|botcoupl=-1;
    expression nor=(nsign*normal(slip));
    nor.write(regionunion({axis,wall,botcoupl}),savefolder+"fixed_normal.vtu");


    //Unknown fields defintion
    //--------------------------

    field u("h1xyz"),v("h1xyz"), p("h1"),T("h1"), y("y"),x ("x"),L("h1xyz");
    //Setting interpolation order :
    u.setorder(membrane, 2);
    v.setorder(KviT, 2);
    p.setorder(fluid, 1);
    T.setorder(KviT, 2);
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
    chambre += integral(KviT,(rho*div(dof(v))+rho*(beta*(dt(dof(p)))-alpha_p*(dt(dof(T)))))*tf(p) );

    //Navier-Stokes
    //--------------
    chambre += integral(KviT,-rho*dt(dof(v))*tf(v));
    chambre += integral(KviT,-grad(dof(p))*tf(v)); 
    chambre += integral(KviT,-lmb*(frobeniusproduct(grad(dof(v)),grad(tf(v)))+frobeniusproduct(transpose(grad(dof(v))),grad(tf(v)))));
    chambre += integral(KviT,(2*lmb/3-muB)*div(dof(v))*trace(grad(tf(v))) );
    //non linearity terms of NS
    //expression(v).reuseit();
    //chambre += integral(fluid,-rho*( grad(v)*dof(v) + grad(dof(v))*v - grad(v)*v )*tf(v));

    //thermal process
    //----------------
    chambre+=integral(KviT,kappa*grad(dof(T))*grad(tf(T)));
    chambre+=integral(KviT,(-alpha_p*T0*dt(dof(p))+rho*Cp*dt(dof(T)))*tf(T));


    //laser heat source in the cavity volume
    //-----------------------------------------
    double shift=.5e-6,pulsewidth=1e-7,waist=100e-6,absp=1,power=1;
    //spatial shape
    expression laser_shape=power*2/(getpi()*pow(waist,2))*pow(2.71828,-2*pow(x/waist,2));
    //gausian time pulse
    expression time_shape=pow(2.71828,-pow(getpi()/pulsewidth*(t()-shift),2));//
    //square time pulse
    //expression time_shape(t()-(shift-pulsewidth/2),expression(-t()+(shift+pulsewidth/2),power/pulsewidth,0),0);
    //Beert-Lambert law for absorption
    chambre+=integral(KviT,-absp*pow(2.71828,absp*y)*laser_shape*time_shape*tf(T));


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
    vec init(chambre);
    //init.setdata(fluid,T);
    //init.setdata(fluid,p);


    // Define the Gen alpha object to time-solve formulation 'chambre' 
    //------------------------------------------------------------------

    genalpha gena(chambre, init, vec(chambre), vec(chambre), {false,true,true,true});
    gena.setparameter(0.75);

    gena.settolerance(1e-10);

    //variable time stepping resolution
    //--------------------------------------
    //create as many vectors as different time steps
    //each time the last timestep is not computed, start the next one with the previous last one
    std::vector<std::vector<double>> timesteping;
    timesteping={

        //{0,10e-9,4e-6,10}
        // different time steps
        {0,.01e-6,4.999e-6,10}
        ,{5e-6,.05e-6,14.99e-6,2}
        ,{15e-6,.1e-6,24.99e-6,1}
        ,{25e-6,.25e-6,99.99e-6,1}
        //,{100e-6,1e-6,1e-3,10}
	
    };


    //parsing for file handling
    //------------------------------------
    //evaluate the total amount of steps to be computed and written
    int tot_steps=0;
    int tot_steps_rec=0;
    double temp;
    for(int i=0; i<timesteping.size() ; i++)
    {
        temp=(timesteping[i][2]-timesteping[i][0])/(timesteping[i][1]);
        tot_steps += 1+floor(temp);
        tot_steps_rec += 1+floor(temp/timesteping[i][3]);
    }
    std::cout<<"total time steps to be computed: "<< std::to_string(tot_steps)<< std::endl;
    std::cout<<"total time steps to be recorded: "<< std::to_string(tot_steps_rec)<< std::endl;

    //create a string with the same number of zeros as digits in the total amount of steps
    std::string init_num_file=std::to_string(tot_steps_rec);
    for(int i=0; i < init_num_file.size() ; i++)
    {
        init_num_file[i]='0';
    }

    //create probe object
    //-----------------------
    //create a record file for the probe and set it

    std::string filename=savefolder+"my_probes.txt";
    ofstream probes;
    probes.open(filename);
    //the 2 first column are reserved
    probes << "step	"<<"time	";
    //define here the names of each of your probes, and use a tab after each name
    probes <<"max p	" <<"max T	"<<"max v_f	"<<"max u_mb	"<<"max v_mb	"<<"norm v	"<<"norm p	"<<"norm T";
    probes << endl;
    probes.close();

    //Create time stepping pvd file which register each vtu file to its own time
    //-----------------------------------------------------------------------------
    std::string filename2=savefolder+"times.pvd";
    ofstream timeorder;
    timeorder.open(filename2);
    //the 2 first column are reserved
    timeorder << "<?xml version=\"1.0\"?>\n";
    timeorder << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    timeorder << "	<Collection>\n";
    probes.close();


    //initiate solution vectors
    std::vector<std::vector<vec>> ssol(timesteping.size());
    std::vector<std::vector<vec>> solution(timesteping.size());//this one is for the field F
    std::vector<std::vector<vec>> solution2(timesteping.size());//that one for dF/dt

    wallclock clk2;
    int index=0;
    //initiate a vector for tracking all the saved time steps
    std::vector<double> times(tot_steps);

    //extra field for display v_membrane
    field dtu("h1xyz");
    dtu.setorder(membrane, 2);


    for (int i=0;i<timesteping.size();i++)
    {

        clk2.tic();
        //Linear or Non-linear solvers
        ssol=gena.runlinear(timesteping[i][0], timesteping[i][1], timesteping[i][2], timesteping[i][3]);
        //ssol=gena.runnonlinear(timesteping[i][0], timesteping[i][1], timesteping[i][2],-1, timesteping[i][3],1);

        //assigning from vector solution for:
        //field F [0], [1] dF/dt   [2] d2F/dt2
        solution[i]=ssol[0];
        solution2[i]=ssol[1];

        std::cout << "File_batch/probes written :"<< std::flush;

        for (int ts = 0; ts < solution[i].size(); ts++)
        {
            // Transfer the data to the field:
            u.setdata(membrane, solution[i][ts]);
            dtu.setdata(membrane,(solution2[i][ts]|u));
            v.setdata(KviT, solution[i][ts]);
            p.setdata(fluid, solution[i][ts]);
            T.setdata(KviT, solution[i][ts]);
            //L.setdata(mechacoucoupl, solution[i][ts]);

            //parsing file number
            std::string indecs=std::to_string(index);
            std::string num_file=init_num_file;
            //file number written
            for (int i=1;i<=indecs.size();i++)
            {
                num_file[num_file.size()-i]=indecs[indecs.size()-i];
            }


            // Write with an order 2 interpolation and with the name of your choice (the name is also the relative path)
            u.write(membrane, savefolder+"u"+num_file+".vtu",2); 
            //dtu.write(membrane, savefolder+"dtu"+num_file+".vtu",2); 
            v.write(KviT, savefolder+"v"+num_file+".vtu",2);		
            p.write(KviT, savefolder+"pch"+num_file+".vtu",2);
            p.write(outair, savefolder+"pout"+num_file+".vtu",2);
            T.write(KviT, savefolder+"T"+num_file+".vtu",2); 

            times[index]=timesteping[i][0]+ts*timesteping[i][3]*timesteping[i][1];

            timeorder.open(filename2, std::ios_base::app);
            timeorder << "		<DataSet timestep=\"" << times[index] << "\" part=\"0\"\n";
            timeorder << "			file=\"pch"+num_file+".vtu\"/>\n";
            timeorder << "		<DataSet timestep=\"" << times[index] << "\" part=\"0\"\n";
            timeorder << "			file=\"pout"+num_file+".vtu\"/>\n";
            timeorder << "		<DataSet timestep=\"" << times[index] << "\" part=\"0\"\n";
            timeorder << "			file=\"T"+num_file+".vtu\"/>\n";
            timeorder << "		<DataSet timestep=\"" << times[index] << "\" part=\"0\"\n";
            timeorder << "			file=\"v"+num_file+".vtu\"/>\n";
            timeorder << "		<DataSet timestep=\"" << times[index] << "\" part=\"0\"\n";
            timeorder << "			file=\"u"+num_file+".vtu\"/>\n";
            timeorder << "		<DataSet timestep=\"" << times[index] << "\" part=\"0\"\n";
            timeorder << "			file=\"dtu"+num_file+".vtu\"/>\n";
            timeorder.close();

            //compute probes values and record it in the file
            //open the file
            probes.open(filename, std::ios_base::app);
            probes.precision(15);
            //step and time index are first
            probes << index <<"  "<< times[index] << "	";
            // put then each expression with the same order than their definition
            probes << expression(p).max(KviT,5)[0] <<"	" ;
            probes << expression(T).max(KviT,5)[0] <<"	" ;
            probes << expression(sqrt(v.compy()*v.compy()+v.compx()*v.compx())).max(KviT,5)[0] <<"	" ;
            probes << norm(v).integrate(KviT, 5) << "	" ;
            probes << norm(p).integrate(KviT, 5) << "	" ;
            probes << norm(T).integrate(KviT, 5) ;
            //probes << expression(sqrt(u.compy()*u.compy()+u.compx()*u.compx())).max(membrane,5)[0]<<"	" ;
            //probes << expression(dtu.compy()).max(membrane,5)[0] ; 
            //end the line
            probes << endl;
            //close the file
            probes.close();

            index++;	

            std::cout <<"  "<< index << std::flush;

        }	
		
		
        std::cout << std::endl;

        clk2.print("step took: ");


    }
    timeorder.open(filename2, std::ios_base::app);
    timeorder << "</Collection>\n";
    timeorder << "</VTKFile>\n";
    timeorder.close();

    clk.print("Total time elapsed:  ");
    probes.open(filename, std::ios_base::app);
    probes << endl << endl <<"Total time elapsed:  "<<  clk.toc()/1e9 << " sec" << endl;
    probes.close();

}

// THE MESH BELOW IS FULLY STRUCTURED AND IS CREATED USING THE (BASIC) SPARSELIZARD GEOMETRY CREATION TOOL.
// BOUNDARY LAYER IS ALSO INCLUDED AS WELL AS A CLOSE LINE TO THE Z-AXIS TO ACT AS THE REVOLUTION AXE

mesh createmesh(double rmb, double depth, double thick, double BLth)
{
    //closeness of the axis
    double caxis=1e-6, rarc=2e-3;
    // Give names to the physical region numbers:
    int membrane = 1, KviT = 2, topair = 3, wall=4, laser=5, BL=6 , clamp=8,axis=9 , outair=10,outrad=11,botcoupl=12, P1=13,P2=14, topcoupl=15; ; 

    // Number of mesh layers:
    int nxmb = 49, nzmb = 4, nzKviT = 73, nair = 50, nBLth=10;

    // define a vector of nodes along a vertical line

    std::vector<double>KviT_vpoints(3*(nzKviT+2*nBLth+1));

    double step1=sqrt(BLth)/nBLth, step2=(depth-2*BLth)/nzKviT;

    // nodes are spaced accordignly to have a boundary layer mesh

    for( int i=0;i<KviT_vpoints.size()/3;i++)
    {
        if(i<=nBLth)
        {
            // if linear BL :/*i*BLth/nBLth*/
            KviT_vpoints[3*i]=0;KviT_vpoints[3*i+1]=-std::pow(i*step1,2);KviT_vpoints[3*i+2]=0;	
            //std::cout<<i<<" "<<-std::pow(i*step1,2)<<std::endl;
        }
        else
        {
            if(i> nBLth && i<=nBLth+nzKviT)
            {
                KviT_vpoints[3*i]=0;KviT_vpoints[3*i+1]=-BLth-(i-nBLth)*step2;KviT_vpoints[3*i+2]=0;	
            }
            else
            {
                // if linear BL :-(i-nzKviT-2*nBLth)*BLth/nBLth
                KviT_vpoints[3*i]=0;KviT_vpoints[3*i+1]=-depth+std::pow((i-nzKviT-2*nBLth)*step1,2);KviT_vpoints[3*i+2]=0;		
            }
        }

        }
    KviT_vpoints[3*(nzKviT+2*nBLth)]=0;KviT_vpoints[3*(nzKviT+2*nBLth)+1]=-depth;KviT_vpoints[3*(nzKviT+2*nBLth)+2]=0;	

    // define a vector of nodes along a horizontal line

    std::vector<double>KviT_hpoints(3*(2+nxmb+nBLth+1));

    step1=sqrt(BLth)/nBLth;
    step2=(rmb-BLth-caxis)/(nxmb);

    //special node for the axis ; should be smaller than one node

    KviT_hpoints[0]=0;KviT_hpoints[1]=0;KviT_hpoints[2]=0;
    KviT_hpoints[3]=caxis;KviT_hpoints[4]=0;KviT_hpoints[5]=0;

    for( int i=2;i<KviT_hpoints.size()/3;i++)
    {
        if(i<=nxmb)
        {
            KviT_hpoints[3*i]=caxis+(i-1)*step2;KviT_hpoints[3*i+1]=0;KviT_hpoints[3*i+2]=0;			
        }
        else
        {
            // if linear BL : (i-2-nxmb-nBLth)*BLth/nBLth
            KviT_hpoints[3*i]=rmb-std::pow((i-2-nxmb-nBLth)*step1,2);KviT_hpoints[3*i+1]=0;KviT_hpoints[3*i+2]=0;	
        }

    }

    // define lines with nodes predefined in previous vectors

    shape ligne1("line",axis,KviT_vpoints);
    shape ligne2("line",botcoupl,KviT_hpoints);
    shape ligne3("line",-1,KviT_vpoints);
    shape ligne4("line",-1,KviT_hpoints);

    ligne3.shift(rmb,0,0);
    ligne4.shift(0,-depth,0);
    //printvector(ligne1.getcoords());
    //printvector(ligne2.getcoords());

    // derive rectangle with boundary layer

    shape nK("quadrangle",KviT,{ligne1,ligne4,ligne3,ligne2});
    //shape nK("quadrangle",KviT,{0,0,0,0,-depth,0,rmb,-depth,0,rmb,0,0},{nzKviT,nxmb,nzKviT,nxmb});

    shape ligne5("line",topcoupl,KviT_hpoints);
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

    //std::cout<<KviT_hpoints.size()/3<<"   "<< 3+nxmb+nBLth<<"    "<<sin(3.14/4)<<endl;

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

