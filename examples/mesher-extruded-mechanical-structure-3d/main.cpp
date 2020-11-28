// This code creates a structured mesh for a 3D geometry of a mechanical part (a small and a bigger cylinder connected together).
//
// It is structured as follows:
//
// 1. The 2D projected surface of the geometry is created
// 2. The projection is extruded
// 3. A deformation function is applied to the big cylinder to get the final geometry.
// 4. A mechanical simulation is performed on the geometry
//
// In order to identify the regions needed in the simulation the following physical region numbers are used:
//
// --> 1 for the whole volume
// --> 2 for the bottom surface of the small cylinder (mechanical clamp)
// --> 3 for the surface inside the big cylinder
// 
// Region number -1 is used for all construction geometries (not used in the simulation).


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The physical region numbers as detailed above:
    int vol = 1, clamp = 2, faceincylinder = 3;
    
    // Define the x, y and z coordinate fields:
    field x("x"), y("y"), z("z");
    
    // Define pi:
    double pi = getpi();
    
    
    ///// CONSTRUCTION OF THE SMALL CYLINDER with inner and outer radius 7mm and 15mm:
    
    double smallinnerradius = 0.007, smallouterradius = 0.015;
    
    // Define the points to construct the cylinder.
    // Inputs are the x, y and z point coordinates.
    shape p1("point", -1, {0,0,0});
    shape p2("point", -1, {smallinnerradius*cos(45*2*pi/360), smallinnerradius*sin(45*2*pi/360), 0});
    shape p3("point", -1, {smallinnerradius*cos(-45*2*pi/360),smallinnerradius*sin(-45*2*pi/360),0});
    shape p4("point", -1, {smallouterradius*cos(45*2*pi/360), smallouterradius*sin(45*2*pi/360), 0});
    shape p5("point", -1, {smallouterradius*cos(-45*2*pi/360),smallouterradius*sin(-45*2*pi/360),0});
    
    // Define the arcs and the lines to construct the small cylinder.
    // Inputs are the first point, last point and arc center point 
    // as well as the number of mesh nodes in the arc.
    shape a1("arc", -1, {p2,p3,p1}, 16);
    shape a2("arc", -1, {p4,p5,p1}, 16);
    shape a3("arc", -1, {p3,p2,p1}, 8);
    shape a4("arc", -1, {p5,p4,p1}, 8);
    
    shape l1("line", -1, {p2,p4}, 3);
    shape l2("line", -1, {p3,p5}, 3);
    
    // Define the 2 quadrangle surfaces for the small cylinder.
    // Inputs are the contour lines (following the contour clockwise or counter-clockwise).
    shape q1("quadrangle", -1, {a1,l1,a2,l2});
    shape q2("quadrangle", -1, {a3,l2,a4,l1});
    
    
    ///// CONSTRUCTION OF THE BIG CYLINDER with inner and outer radius 10mm and 30mm
    // around center point on x axis at x equal 100mm.
    
    double biginnerradius = 0.01, bigouterradius = 0.03;
    
    // Define the points to construct the cylinder.
    shape p6("point", -1, {0.1,0,0});
    shape p7("point", -1, {0.1-biginnerradius*cos(45*2*pi/360), biginnerradius*sin(45*2*pi/360), 0});
    shape p8("point", -1, {0.1-biginnerradius*cos(-45*2*pi/360),biginnerradius*sin(-45*2*pi/360),0});
    shape p9("point", -1, {0.1-bigouterradius*cos(45*2*pi/360), bigouterradius*sin(45*2*pi/360), 0});
    shape p10("point", -1, {0.1-bigouterradius*cos(-45*2*pi/360),bigouterradius*sin(-45*2*pi/360),0});
    
    // Define the arcs and the lines to construct the big cylinder.
    shape a5("arc", -1, {p7,p8,p6}, 8);
    shape a6("arc", -1, {p9,p10,p6}, 8);
    shape a7("arc", -1, {p8,p7,p6}, 16);
    shape a8("arc", -1, {p10,p9,p6}, 16);
    
    shape l3("line", -1, {p7,p9}, 5);
    shape l4("line", -1, {p8,p10}, 5);
    
    // Define the 2 quadrangle surfaces for the big cylinder.
    shape q3("quadrangle", -1, {a5,l3,a6,l4});
    shape q4("quadrangle", -1, {a7,l4,a8,l3});
    
    
    ///// CONSTRUCTION OF THE SURFACE BETWEEN THE TWO CYLINDERS
    
    // Define the lines linking the 2 cylinders:
    shape l5("line", -1, {p4,p9}, 12);
    shape l6("line", -1, {p5,p10}, 12);
    
    // Define the surface between the cylinders:
    shape q5("quadrangle", -1, {a4,l6,a6,l5});
    
    
    ///// EXTRUDE THE TOP AND BOTTOM CYLINDERS AROUND THE CENTRAL PART OF THE SMALL CYLINDER:
    // Extrude arguments are the physical region of the extruded volume, 
    // the extrusion height and the numbers of node layers in the extrusion.
    
    // The central part in the small cylinder is 10mm thick and the two above and below 5mm.
    double centralthickness = 0.01, thicknessaboveandbelow = 0.005;
    
    // Create the bottom volume part:
    shape v1 = q1.extrude(vol, -thicknessaboveandbelow, 3);
    shape v2 = q2.extrude(vol, -thicknessaboveandbelow, 3);
    
    // Create the central part:
    shape v3 = q1.extrude(vol, centralthickness, 6);
    shape v4 = q2.extrude(vol, centralthickness, 6);
    
    // Create the top volume part.
    // Extrude it from the top face of v3 and v4.
    shape v3topface = v3.getsons()[5];
    shape v4topface = v4.getsons()[5];
    
    shape v5 = v3topface.extrude(vol, thicknessaboveandbelow, 3);
    shape v6 = v4topface.extrude(vol, thicknessaboveandbelow, 3);
    
    
    ///// EXTRUDE THE BIG CYLINDER:
    
    shape v7 = q3.extrude(vol, centralthickness, 6);
    shape v8 = q4.extrude(vol, centralthickness, 6);
    
    
    ///// EXTRUDE THE SURFACE BETWEEN THE 2 CYLINDERS:
    
    shape v9 = q5.extrude(vol, centralthickness, 6);
    
    
    ///// DEFORM THE BIG CYLINDER ABOVE AND BELOW:
    
    v7.move(array3x1(0,0, 100*(bigouterradius-sqrt((x-0.1)*(x-0.1)+y*y))*(z-centralthickness/2)));
    v8.move(array3x1(0,0, 100*(bigouterradius-sqrt((x-0.1)*(x-0.1)+y*y))*(z-centralthickness/2)));
    
    
    ///// GET THE BOTTOM SURFACE OF THE SMALL CYLINDER AND ASSIGN PHYSICAL REGION NUMBER 2
    ///// GET THE SURFACE INSIDE THE BIG CYLINDER AND ASSIGN PHYSICAL REGION NUMBER 3
    
    // You can .write the mesh at any time during the geometry construction 
    // phase (see below) to easily know which son you are looking for.
    shape clampface1 = v1.getsons()[5];
    shape clampface2 = v2.getsons()[5];
    
    clampface1.setphysicalregion(clamp);
    clampface2.setphysicalregion(clamp);
    
    shape faceincylinder1 = v7.getsons()[1];
    shape faceincylinder2 = v8.getsons()[1];
    
    faceincylinder1.setphysicalregion(faceincylinder);
    faceincylinder2.setphysicalregion(faceincylinder);
    
    
    ///// LOAD THE MESH
    
    // Add here all regions needed in the finite element simulation.
    mesh mymesh({v1,v2,v3,v4,v5,v6,v7,v8,v9, clampface1,clampface2, faceincylinder1,faceincylinder2});
    
    // You can write the mesh at any time during the geometry 
    // construction to easily debug and validate every line of code.
    mymesh.write("mymesh.msh");
    
    
    ///// START THE MECHANICAL SIMULATION
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the mechanical displacement.
    field u("h1xyz");
    
    // Use interpolation order 2 on 'vol', the whole domain:
    u.setorder(vol, 2);
    
    // Clamp on surface 'clamp' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);
    
    // E is Young's modulus. nu is Poisson's ratio. 
    double E = 200e9, nu = 0.3;
    
    formulation elasticity;
    
    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), E, nu));
    // Add a surface force (to the inside of the big cylinder) to tilt the structure around the x axis:
    elasticity += integral(faceincylinder, array1x3(0,10*(z-centralthickness/2),0)*tf(u));
    
    elasticity.generate();
    
    vec solu = solve(elasticity.A(), elasticity.b());
    
    // Transfer the data from the solution vector to the u field:
    u.setdata(vol, solu);
    // Write the deflection with an order 2 interpolation. Exaggerate the deflection by a large factor.
    (0.3e10*u).write(vol, "u.pos", 2);

    // Code validation line. Can be removed.
    std::cout << (abs(compz(u)).max(vol, 2)[0] < 4.24272e-12 && abs(compz(u)).max(vol, 2)[0] > 4.2427e-12);
}

