// This code shows how to perfom a magnetostatic analysis on a 2D cross-section of a
// rotating PMSM (permanent magnet synchronous) electric motor. The motor torque is
// calculated based on the virtual work principle for an increasing mechanical angle.
//
// Antiperiodicity is used to reduce the computational domain to only 45 degrees of
// the total geometry (the motor has 4 pairs of poles). In order to link the rotor and
// stator domains at their interface a general mortar-based continuity condition is used.
// This allows to work with the non-matching mesh at the interface when the rotor moves.

#include <fstream>
#include <experimental/filesystem>
#include <iostream>
#include "sparselizard.h"
#include <vector>


#define WRITE_B false

using namespace std;
using namespace sl;

expression normalize(expression);
expression polar(expression);
expression polarscalar(expression);

expression projection(expression e, expression a);

const char* getVar(int argc, char* argv[], int position) {
    if(argc > position) {
        return argv[position];
    } else {
        cout
            << "Use this software as follows:" << endl
        << "Example command to load only one mesh file" << endl << endl
        << "./loadem.sh FilenameToLoad" << endl << endl
        << "Example command to load only one mesh file save to folderName (default:December)" << endl << endl
        << "./loadem.sh FilenameToLoad FolderName" << endl << endl
        << "Example command to load only multiple mesh file of the name convention: FilenameToLoad{alpha}A{ironAlpha}.msh" << endl << endl << endl
        << "./loadem.sh FilenameToLoad alphaStep ironStep pieces ironpieces" << endl << endl;
        exit(1);
    }
}
void writeVec(vector<double> v, string filename) {
    ofstream myfile(filename);
    int vsize = v.size();
    for (int n=0; n<vsize; n++)
    {
        myfile << v[n] << endl;
    }
}

string twodigits(int num) {
    return ((num <= 9) ? "0" : "") + (to_string((int)(num)));
}



void createPath(string path) {
    string command = "mkdir -p " + path;
    bool success = system(command.c_str());
    // namespace fs = experimental::filesystem; // In C++17 use filesystem.

    // error_code ec;
    // bool success = fs::create_directories(path, ec);

    if (!success) {
        //cout << "Couldn't create directory path: " << path << endl; // Fun fact: In case of success ec.message() returns "The operation completed successfully." using vc++.
    }
}

void writeFile(string filename, string data) {
    ofstream myfile;
    myfile.open (filename.c_str());
    myfile << data;
    myfile.close();
}

struct DataArray {
    vector<double> array;
    string name;

    DataArray(string n) {

        this->name = n;
    }

    void add(double x) {
        this->array.push_back(x);
    }

    void save(string filename) {
        string path = "results/" + this->name;
        createPath(path);
        createPath(path + "/Arrays");
        writeVec(this->array, path + "/Arrays/" + filename + ".csv");
    }
};

struct DataPack {
    expression A;
    expression B;
    expression H;
    expression F;
    double T;
    double TBottom;
    double TTop;
    double TIron;
    double TBM;
    double TTM;
    double alpha;
    double ironAlpha;

    string name;
    int allRegion;

    DataPack(string name) {
        this->name = name;
    }

    void saveData() {
        cout << "Saving datapack to folder " << this->name << endl;

        string path = "results/" + this->name;
        createPath(path);

        createPath(path + "/B");
        createPath(path + "/A");
        createPath(path + "/F");
        createPath(path + "/T");
        createPath(path + "/H");
        createPath(path + "/Iterations");

        int all = this->allRegion;
        int alpha = this->alpha;
        int ironAlpha = this->ironAlpha;

//        this->A.write(all, path + "/A/AField"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        if(WRITE_B) {
            this->B.write(all, path + "/B/BField"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        }
//        this->H.write(all, path + "/H/HField"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
//        this->F.write(all, path + "/F/Force"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        string data = "Bottom,Iron,Top,Bottom Magnets,Top Magnets\r\n" +
            to_string(this->TBottom) + "," + to_string(this->TIron) + "," + to_string(this->TTop) + "," + to_string(this->TBM) + "," + to_string(this->TTM);
        writeFile(path + "/T/torgues"+twodigits(alpha)+twodigits(ironAlpha)+".csv", data);
        // polar(this->B).write(all, path + "TestPolarB"+twodigits(alpha)+".vtu", 2);
        // polar(this->F).write(all, path + "TestPolarF"+twodigits(alpha)+".vtu", 2);
        // polarscalar(this->B).write(all, path + "TestPolarScalarB"+twodigits(alpha)+".vtu", 2);
    }
};


// Input is rotor angular position in degrees. Output is torque in Nm.
DataPack* sparselizard(mesh mymesh, double alpha = 0, double ironAlpha = 0, bool autoload = false, string name="December", bool timeit = false)
{
    wallclock clk;
    // Magnet Power
    double MagnetB = 1.3; // N45 = [1.32-1.38]

    int bottomStator = 1,
        firstmagnetup = 2,
        firstmagnetdown = 3,
        ironbars = 4,
        secondmagnetup = 6,
        secondmagnetdown = 7,
        air = 10,
        topStator = 20,
        bottomSurface = 201,
        topSurface = 202,
        nonmagnetic = 300,
        environment = 301,
        boundary = 400;


    // Define new physical regions for convenience:
    int all = selectall();
    int magnetsdown = selectunion({secondmagnetdown, firstmagnetdown});
    int magnetsup = selectunion({secondmagnetup, firstmagnetup});
    int magnets = selectunion({firstmagnetdown, firstmagnetup, secondmagnetdown, secondmagnetup});
    int iron = selectunion({ironbars, bottomStator, topStator});

    int noniron = selectunion({ nonmagnetic, magnets });
    int magnetic = selectunion({magnets, iron});
    // parts
    int topPart = selectunion({ secondmagnetup, secondmagnetdown, topStator });
    int bottomPart = selectunion({ firstmagnetup, firstmagnetdown, bottomStator });
    int rotor = selectunion({secondmagnetup, secondmagnetdown});
    // sections
    int bottommagnets = selectunion({firstmagnetup, firstmagnetdown});
    int topmagnets = selectunion({secondmagnetup, secondmagnetdown});
    int middlesection = selectunion({ironbars});

    cout << "Calculating Spanning Tree" << endl;
    // Define a spanning tree to gauge the magnetic vector potential (otherwise the matrix is singular).
    // Start growing the tree from the regions with constrained potential vector (here the contour):
    spanningtree spantree({ boundary });

    field A("hcurl", spantree),
         x("x"), y("y"), z("z");
    // gauge the macnetic vector potential
    A.setgauge(all);
    A.setorder(all, 2); // Use interpolation order 2:
    // Put a magnetic wall
    //A.setconstraint(topSurface);

    expression B = curl(A);
    // The remanent induction field in the magnet is 0.5 Tesla perpendicular to the magnet:
    expression normedradialdirection = array3x1(0,0,1);
    expression bremanent = MagnetB * normedradialdirection;
    expression minusnormedradialdirection = array3x1(0,0,-1);
    expression minusbremanent = MagnetB * minusnormedradialdirection;

    // Define the permeability in all regions.
    // Taking into account saturation and measured B-H curves can be easily done
    // by defining an expression based on a 'spline' object (see documentation).
    // BH Curve first try
//    std :: vector<double> Bcurve = {0.0, 0.2,   0.426, 0.761,  1.097, 1.233,  1.335,   1.460,   1.590,   1.690,   1.724,   1.740 , 2.0, 2.5, 10.0, 100.0};
//    std :: vector<double> Hcurve = {0.0, 318.3, 477.5, 795.8, 1591.6, 2387.3, 3978.9, 7957.8, 15915.5, 31831.0, 44456.3, 55704.3 , 58000.0, 60000.0, 70000.0, 80000.0};
    std::vector<double> Hcurve = {
    0.0000e+00, 5.5023e+00, 1.1018e+01, 1.6562e+01, 2.2149e+01, 2.7798e+01, 3.3528e+01,
    3.9363e+01, 4.5335e+01, 5.1479e+01, 5.7842e+01, 6.4481e+01, 7.1470e+01, 7.8906e+01,
    8.6910e+01, 9.5644e+01, 1.0532e+02, 1.1620e+02, 1.2868e+02, 1.4322e+02, 1.6050e+02,
    1.8139e+02, 2.0711e+02, 2.3932e+02, 2.8028e+02, 3.3314e+02, 4.0231e+02, 4.9395e+02,
    6.1678e+02, 7.8320e+02, 1.0110e+03, 1.3257e+03, 1.7645e+03, 2.3819e+03, 3.2578e+03,
    4.5110e+03, 6.3187e+03, 8.9478e+03, 1.2802e+04, 1.8500e+04, 2.6989e+04, 3.9739e+04,
    5.9047e+04, 8.8520e+04, 1.3388e+05, 2.0425e+05, 3.1434e+05, 4.8796e+05, 7.6403e+05};
    std::vector<double> Bcurve = {
    0.0000e+00, 5.0000e-02, 1.0000e-01, 1.5000e-01, 2.0000e-01, 2.5000e-01, 3.0000e-01,
    3.5000e-01, 4.0000e-01, 4.5000e-01, 5.0000e-01, 5.5000e-01, 6.0000e-01, 6.5000e-01,
    7.0000e-01, 7.5000e-01, 8.0000e-01, 8.5000e-01, 9.0000e-01, 9.5000e-01, 1.0000e+00,
    1.0500e+00, 1.1000e+00, 1.1500e+00, 1.2000e+00, 1.2500e+00, 1.3000e+00, 1.3500e+00,
    1.4000e+00, 1.4500e+00, 1.5000e+00, 1.5500e+00, 1.6000e+00, 1.6500e+00, 1.7000e+00,
    1.7500e+00, 1.8000e+00, 1.8500e+00, 1.9000e+00, 1.9500e+00, 2.0000e+00, 2.0500e+00,
    2.1000e+00, 2.1500e+00, 2.2000e+00, 2.2500e+00, 2.3000e+00, 2.3500e+00, 2.4000e+00};
    spline hbcurve(Bcurve, Hcurve);
    // Define nu as h = nu * b and fix the 0/0 division for the first entry:
    std::vector<double> nudata(Hcurve.size());
    for (int i = 1; i < Hcurve.size(); i++) {
        nudata[i] = Hcurve[i]/Bcurve[i];
    }
    nudata[0] = nudata[1];
    // Use a cubic spline interpolation for nu(b):
    spline nucurve(Bcurve, nudata);

    // Vacuum magnetic permeability [H/m]:
    double mu0 = 4.0*getpi()*1e-7;
    parameter mu;
    expression BHCurve(nucurve, norm(B));
    expression mIron = 1/BHCurve;
    mu|all = mu0;
    mu|magnets = mu0;
    mu|iron = 100*mu0; // mIron

    // dhdb expression based on H-B spline
    expression H = B / mu;
    expression dhdb(hbcurve.getderivative(), norm(B));

    cout << "Formulation" << endl;
    formulation magnetostatics;
    // The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
    magnetostatics += integral(all, 1/mu* curl(dof(A)) * curl(tf(A)) );

    // Steel - this term is used for iron saturation
    // magnetostatics += integral(iron, (H + dhdb * (curl(dof(A)) - B)) * curl(tf(A)));
    // Magnets - the magnetization of the magnets (up and down)
    magnetostatics += integral(magnetsdown, -1/mu* minusbremanent * curl(tf(A)));
    magnetostatics += integral(magnetsup, -1/mu* bremanent * curl(tf(A)));

    int numdofs = magnetostatics.countdofs();
    cout << "Unknowns in formulation: " << endl;
    cout << numdofs << endl;

    if(timeit) clk.print("Setup experiment:\t");
    if(autoload == true) {
        cout << "Autoloading cache..." << endl;
        A.loadraw("cache/aaz"+twodigits(alpha)+twodigits(ironAlpha)+".slz.gz", true);
        cout << "Cache loaded" << endl;
        if(timeit) clk.print("CacheLoaded:\t");
    } else {
        cout << "Solving..." << endl;

/*
///////////////////////////////////////////////
// NON LINEAR SOLUTION
///////////////////////////////////////////////
    // Initial solution is a = 0 (thus b = 0):
    vec x(magnetostatics);

    // Nonlinear iteration:
    double relres = 1, maxb; int iter = 0;
    // disable
    while (relres > 1e-5)
    {
        magnetostatics.generate();

        // Get the A and rhs in A*x = rhs:
        mat Amat = magnetostatics.A();
        vec rhs = magnetostatics.b();

        // Calculate the relative residual:
        relres = (rhs - Amat*x).norm()/rhs.norm();
        std::cout << "Relative residual @" << iter << " is " << relres;

        // Relaxation = including progressively the new solution using:
        // xnew = a*xold + (1-a)*solve(A, b)
        double relaxation = 0.;
        x = solve(Amat, rhs);

        // Update the field a with the solution x:
        setdata(x);

        // Do not allow b to go out of the provided [0, 2.4] T data range for the h(b) curve:
        maxb = norm(B).max(all, 5)[0];
        std::cout << " (max b is " << maxb << " T)" << std::endl;
        if (maxb > 2.4)
            x = 2.4/maxb * x;

        setdata(x);

        iter++;
    }
*/
        // Linear Solution
        magnetostatics.solve();
        //az.writeraw(all, "cache/aaz"+twodigits(alpha)+twodigits(ironAlpha)+".slz.gz", true);
        cout << "Solved!" << endl;
        if(timeit) clk.print("Solving time:\t");
    }

    field magforce("h1xyz");
    magforce.setorder(all, 2);

    // The magnetic force is projected on field 'magforce' on the solid stator region.
    // This is done with a formulation of the type dof*tf - force calculation = 0.
    formulation forceprojection;
    forceprojection += integral(all, dof(magforce)*tf(magforce));
    forceprojection += integral(all, -predefinedmagnetostaticforce(tf(magforce, all), 1/mu*B, mu));
    forceprojection.solve();

    if(timeit) clk.print("Forces calc time:\t");
    // expression magforcescalar = polarscalar((expression)magforce);
    // magforcescalar.write(bottomPart, "results/bottomMagForce"+twodigits(alpha)+".vtu");

    // Calculate the torque:
    expression leverarm = array3x1(x,y,0);
    expression torgue = compz(crossproduct(leverarm, magforce));

    double T = torgue.integrate(bottomPart, 4);
    double TBottom = T;
    double TTop = torgue.integrate(topPart, 4);
    double TIron = torgue.integrate(middlesection, 4);
    double TBM = torgue.integrate(bottommagnets, 4);
    double TTM = torgue.integrate(topmagnets, 4);

    double scaleFactor = 1.;
    cout << "--------------------------------------------------------" << endl;
    cout << "TorgueBottom: " << TBottom * scaleFactor << endl;
    cout << "TorgueIron: " << TIron * scaleFactor << endl;
    cout << "TorgueTop: " << TTop * scaleFactor << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;

    DataPack* pack = new DataPack(name);
    pack->T = T* scaleFactor;
    pack->TBottom = T* scaleFactor;
    pack->TTop = TTop* scaleFactor;
    pack->TIron = TIron* scaleFactor;
    pack->TBM = TBM* scaleFactor;
    pack->TTM = TTM* scaleFactor;

    pack->F = (expression)magforce;
//    pack->A = (expression)az;
    pack->B = curl(A);
//    pack->H = curl(az)/mu;
    pack->alpha = alpha;
    pack->ironAlpha = ironAlpha;
    pack->allRegion = selectall();

    // pack->saveData();
    if(timeit) clk.print("Pack Ready:\t");
    return pack;
}

string getFilename(int alpha,int ironAlpha,int argc,char* argv[], bool autoname = false) {
    // FILE MESH 2.2 .msh
    string num = twodigits(alpha);
    string numiron = twodigits(ironAlpha);
    // string test = "models/New.msh";
    string filename;
    if ( argc < 4 && autoname == false ) {
        // take file from args
        filename = (string)"models/" + argv[1];
        cout << "Loading model: " << filename << endl;
    } else {
        filename = (string)"models/parts/" + argv[1] + num + "A" + numiron + ".msh";
        cout << "Loading model(autoname): " << filename << endl;
    }
    return filename;
}

void simulate(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Provide the name of the mesh file!" << endl;
        return;
    }

    string filename = getFilename(0, 0, argc, argv, true);
    mesh mymesh(filename);

    DataPack* pack;
    pack = sparselizard(mymesh, 0, 0, false);
    pack->saveData();
}

int main(int argc, char *argv[])
{
    PetscInitialize(0,{},0,0);
    wallclock clk;

    if(argc < 2) {
        cout << "Need to specify a model name. e.g. ./mcvt Parrot" << endl;
        return 0;
    }


    int start = 7; // last 7
    int end = 45;
    int alphaStep = atoi(getVar(argc, argv, 2));
    int pieces = (end-start)/alphaStep;


    string folderName = "October";

    createPath("cache");

    int ironstart = 15; // last 15
    int ironend = 45;
    int ironStep = atoi(getVar(argc, argv, 3));
    int ironpieces = (ironend - ironstart)/ironStep;
    cout << "ROT  Step set to: " << alphaStep << endl;
    cout << "IRON Step set to: " << ironStep << endl;

    if (argc >= 6) {
        start = atoi(argv[4]);
        cout << "START ROT  SET TO : " << start << endl;
        ironstart = atoi(argv[5]);
        cout << "START IRON SET TO : " << ironstart << endl;
    }

    DataPack* pack;
    // make only one loop
    for (double ironAlpha = ironstart; ironAlpha <= ironstart+(ironpieces*ironStep); ironAlpha += ironStep)
    {
        DataArray* torguesB = new DataArray(folderName);
        DataArray* torguesI = new DataArray(folderName);
        DataArray* torguesT = new DataArray(folderName);
        for (double alpha = start; alpha <= start+(pieces*alphaStep); alpha += alphaStep)
        {
            cout << "Iron Angle\t" << ironAlpha << "\t Rotor Angle:\t" << alpha << endl;
            string filename = getFilename(alpha, ironAlpha, argc, argv, true);
            mesh mymesh;
            mymesh.selectskin(400); // 400 = boundary
            mymesh.load(filename);
            pack = sparselizard(mymesh, alpha, ironAlpha, false, folderName, true);
            pack->saveData();
            torguesB->add(pack->TBottom);
            torguesI->add(pack->TIron);
            torguesT->add(pack->TTop);
            delete pack;
        }
        torguesB->save("TBottom" + twodigits(ironAlpha));
        torguesI->save("TIron" + twodigits(ironAlpha));
        torguesT->save("TTop" + twodigits(ironAlpha));
        delete torguesB, torguesI, torguesT;
        start = 0;
    }

    clk.print("Total run time:");
    PetscFinalize();
    return 0;
}

expression normalize(expression e)
{
    return e / norm(e);
}

expression polar(expression e) {
    field x("x"), y("y"), z("z");
    // expression proj = e * array3x1(-y,x,0);
    expression proj = projection(e, array3x1(-y,x,0));
    return proj;
}

expression polarscalar(expression e) {
    field x("x"), y("y"), z("z");
    expression proj = projection(e, array3x1(-y,x,0));
    return norm(proj);
}

expression projection(expression e, expression a) {
    return ((e * a) / (a * a)) * a;
}
