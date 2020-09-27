// The purpose of this example is to check that saving to and loading from .slz works as expected when the software is modified.


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(bool iswritemode)
{	
    if (iswritemode)
    {
        int surleft = 1, surright = 2, left = 3, right = 4;
        
        int n = 5;
        shape t1("triangle", surleft, {0,0,0, 1,0,0, 0,1,0}, {n,n,n});
        shape t2("triangle", surright, {1,0,0, 1,1,0, 0,1,0}, {n,n,n});
        shape l1 = t1.getsons()[2];
        l1.setphysicalregion(left);
        shape l2 = t2.getsons()[0];
        l2.setphysicalregion(right);
        mesh mymesh({t1,t2,l1,l2});
        int all = regionunion({surleft,surright});

        field u("h1xy", {2,3}), x("x"), y("y");
        
        u.setorder(surleft, 4);
        u.setorder(surright, 1);
        
        u.harmonic(2).setvalue(surleft, array2x1(x,2));
        u.harmonic(2).setvalue(surright, array2x1(x*x*x,y*y));
        u.harmonic(3).setvalue(all, array2x1(1,2));
        
        double focheck = fieldorder(u.compx()).integrate(all, 5);
        double ucheck = norm(u.harmonic(2)).integrate(all, 5) + norm(u.harmonic(3)).integrate(all, 5);
        
        std::vector<double> extradata = {focheck, ucheck};
        
        mymesh.write("mesh.msh");
        u.writeraw(all, "u.slz", false, extradata);
    }
    else
    {
        int all = 5;
        
        mesh mymesh("mesh.msh");
        
        field u("h1xy", {2,3});
        
        std::vector<double> extradata = u.loadraw("u.slz");
        
        double focheck = fieldorder(u.compx()).integrate(all, 5);
        double ucheck = norm(u.harmonic(2)).integrate(all, 5) + norm(u.harmonic(3)).integrate(all, 5);
        
        std::cout << focheck << " " << ucheck << std::endl;
        std::cout << std::abs(focheck-extradata[0])/std::abs(focheck) << " " << std::abs(ucheck-extradata[1])/std::abs(ucheck) << std::endl;
        
        std::cout << (std::abs(focheck-extradata[0])/std::abs(focheck) == 0 && std::abs(ucheck-extradata[1])/std::abs(ucheck) == 0);
    }
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard(true);
    sparselizard(false);
    
    SlepcFinalize();

    return 0;
}

