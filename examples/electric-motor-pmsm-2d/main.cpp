// This code shows how to perfom a magnetostatic analysis on a 2D cross-section of a
// rotating PMSM (permanent magnet synchronous) electric motor. The motor torque is 
// calculated based on the virtual work principle for an increasing mechanical angle.
//
// Anti-periodicity is used to reduce the computational domain to only 45 degrees of 
// the total geometry (the motor has 4 pairs of poles). In order to link the rotor and
// stator domain at their interface a general Mortar-based continuity condition is used.
// This allows to work with the non-matching mesh at the interface when the rotor moves.  


#include "sparselizardbase.h"


using namespace mathop;

// Input is rotor angular position in degrees. Output is torque in Nm.

double sparselizard(double alpha)
{	
    // Give names to the physical region numbers :
    int rotmagmat = 1, magnet = 2, magnetgap = 3, gaprot = 4, gapstat = 5, statmagmat = 6, windslot = 7, winda = 8, windb = 9, windc = 10;
    int gamma1 = 11, gamma2 = 12, inarc = 13, outarc = 14;

    // Load the mesh. Set verbosity to 0.
    mesh mymesh("pmsm.msh", 0);

    // Define new physical regions for convenience:
    int rotor = regionunion({rotmagmat, magnet, magnetgap, gaprot});
    int windings = regionunion({winda, windb, windc});
    int stator = regionunion({gapstat, windslot, statmagmat, windings});
    int nonmag = regionunion({magnet,magnetgap, gaprot, gapstat, windslot, windings});
    int rotstatinterface = regionintersection({rotor,stator});

    // Define the number of pole pairs (4 in the geometry used):
    int numpolepairs = 4;

    // Periodicity regions gamma1 and gamm2 for the rotor and the stator:
    int rotgamma1 = regionintersection({rotor, gamma1});
    int rotgamma2 = regionintersection({rotor, gamma2});
    int statgamma1 = regionintersection({stator, gamma1});
    int statgamma2 = regionintersection({stator, gamma2});

    // Peak winding current [A] times number of turns:
    double Imax = 300;
    // Calculate the area of a winding:
    double windarea = expression(1).integrate(winda, 5);


    // Nodal shape functions 'h1' for the z component of the vector potential.
    field azrot("h1"), azstat("h1"), x("x"), y("y"), z("z");

    // In 2D the vector potential only has a z component:
    expression arot = array3x1(0,0,azrot), astat = array3x1(0,0,azstat);

    // Use interpolation order 2:
    azrot.setorder(rotor, 2);
    azstat.setorder(stator, 2);

    // Put a magnetic wall at the inner rotor and outer stator boundaries:
    azrot.setconstraint(inarc);
    azstat.setconstraint(outarc);

    // The remanent induction field in the magnet is 1.0 Tesla perpendicular to the magnet:
    expression normedradialdirection = array3x1(x,y,0)/sqrt(x*x+y*y);
    expression bremanent = 0.5 * normedradialdirection;     


    // Vacuum magnetic permeability [H/m]: 
    double mu0 = 4.0*getpi()*1e-7;

    // Define the permeability in all regions.
    //
    // Taking into account saturation and measured B-H curves can be easily done
    // by defining an expression based on a 'spline' object (see documentation).
    //
    parameter mu;

    mu|rotor = 2000.0*mu0;
    mu|stator = 2000.0*mu0;
    // Overwrite on non-magnetic regions:
    mu|nonmag = mu0;


    formulation magnetostatics;

    // The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
    magnetostatics += integral(rotor, 1/mu* curl(dof(arot)) * curl(tf(arot)) );
    magnetostatics += integral(stator, 1/mu* curl(dof(astat)) * curl(tf(astat)) );

    // Add the remanent magnetization of the rotor magnet:
    magnetostatics += integral(magnet, -1/mu* bremanent * curl(tf(arot)));

    // Add the current density source js [A/m2] in the z direction in the stator windings.
    // A three-phased actuation is used. The currents are dephased by the mechanical angle
    // times the number of pole pairs. This gives a stator field rotating at synchronous speed.

    // Change the phase (degrees) to tune the electric angle: 
    double phase = 0;

    parameter jsz;

    jsz|winda = Imax/windarea * sin( (phase + 4.0*alpha - 0.0) * getpi()/180.0);
    jsz|windb = Imax/windarea * sin( (phase + 4.0*alpha - 60.0) * getpi()/180.0);
    jsz|windc = Imax/windarea * sin( (phase + 4.0*alpha - 120.0) * getpi()/180.0);

    magnetostatics += integral(windings, array3x1(0,0,jsz) * tf(astat));

    // Rotor-stator continuity condition (including antiperiodicity settings with factor '-1'):
    magnetostatics += continuitycondition(rotstatinterface, rotstatinterface, azrot, azstat, {0,0,0}, alpha, 45, -1);

    // Rotor and stator antiperiodicity condition:
    magnetostatics += periodicitycondition(rotgamma1, rotgamma2, azrot, {0,0,0}, {0,0,45}, -1);
    magnetostatics += periodicitycondition(statgamma1, statgamma2, azstat, {0,0,0}, {0,0,45}, -1);


    solve(magnetostatics);


    std::string alphastring = std::to_string((int)alpha);
    azstat.write(stator, "astator"+alphastring+".pos", 2);
    curl(astat).write(stator, "bstator"+alphastring+".pos", 2);

    // Write the rotor fields on the rotated mesh:
    mymesh.rotate(0,0,alpha);

    azrot.write(rotor, "arotor"+alphastring+".pos", 2);
    curl(arot).write(rotor, "brotor"+alphastring+".pos", 2);

    // Rotate back to the original position:
    mymesh.rotate(0,0,-alpha);


    // The MAGNETOSTATIC FORCE acting on the rotor is computed below.

    // This field will hold the x and y component of the magnetic forces:
    field magforce("h1xy");

    // The magnetic force is projected on field 'magforce' on the solid rotor region.
    // This is done with a formulation of the type dof*tf - force calculation = 0.
    formulation forceprojection;

    forceprojection += integral(statmagmat, dof(magforce)*tf(magforce));
    expression Hstat = 1/mu * (curl(astat));
    forceprojection += integral(stator, - predefinedmagnetostaticforce(tf(magforce, statmagmat), Hstat, mu));

    solve(forceprojection);

    // Calculate the torque:
    expression leverarm = array3x1(x,y,0);

    double torque = compz(crossproduct(leverarm, magforce)).integrate(statmagmat, 5);

    // The torque has to be scaled to the actual motor z length (50 mm) and multiplied
    // by 8 to take into account the full 360 degrees of the motor.
    // A minus sign gives the torque on the rotor (opposite sign than the stator torque).
    torque = - torque * 8.0 * 0.05;

    return torque;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    wallclock clk;

    std::cout << "Mechanical angle [degrees] and torque [Nm]:" << std::endl;
    for (double alpha = 0.0; alpha <= 45.0; alpha += 1.0)
    {
        double torque = sparselizard(alpha);
        std::cout << alpha << " " << torque << std::endl;   
    }

    clk.print("Total run time:");

    SlepcFinalize();

    return 0;
}

