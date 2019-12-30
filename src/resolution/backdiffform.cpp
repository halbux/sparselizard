#include "backdiffform.h"

backdiffform::backdiffform(formulation formul, int ordre, vec xinit, vec dtxinit, std::vector<bool> isrhskcconstant)
{
    myformulation = formul;
    
    x = xinit;
    dtx = dtxinit;
	order=ordre;
    
    if (isrhskcconstant.size() == 0)
        isconstant = {false,false,false};
    else
        isconstant = isrhskcconstant;
    if (isconstant.size() != 3)
    {
        std::cout << "Error in 'backdiffform' object: expected a length 3 or empty vector as fourth argument" << std::endl;
        abort();  
    }
	if( order < 1 || order > 6)
	{
        std::cout << "Error in 'backdiffform' order: expected a non-zero integer up to 6. Above, the method is unstable." << std::endl;
        abort();  
    }	
}

void backdiffform::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void backdiffform::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }

std::vector<std::vector<vec>> backdiffform::runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep, int verbosity)
{
    return run(true, starttime, timestep, endtime, -1, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> backdiffform::runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    return run(false, starttime, timestep, endtime, maxnumnlit, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> backdiffform::run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    // Solve end time rounding issues:
    endtime += endtime*1e-12;
    
    if (starttime > endtime)
        return {};
    
    if (outputeverynthtimestep <= 0)
        outputeverynthtimestep = 1;
        
    // Get all fields in the formulation:
    std::shared_ptr<dofmanager> dofmngr = myformulation.getdofmanager();
    std::vector<std::shared_ptr<rawfield>> allfields = dofmngr->getfields();
    // Set all fields in the formulation to the initial solution:
    for (int i = 0; i < allfields.size(); i++)
        allfields[i]->setdata(-1, x|field(allfields[i]));
    
    vec rhs, right_mbr; mat K, C, leftmat;
	
	std::vector<vec> solution_pool(order);
	for(int i=0;i<order;i++)
		solution_pool[i]=x;
    
    // Count the number of time steps to step through and the number of vectors to output:
    int numtimesteps = 0; int outputsize = 0;
    for (double t = starttime; t <= endtime; t = t + timestep)
    {
        if (numtimesteps%outputeverynthtimestep == 0)
            outputsize++;
        numtimesteps++;
    }
    
    
    // Start the Backward Differetiation Formula scheme iteration:
    std::cout << "Backward Differetiation Formula with order :"<<order<<"   for " << numtimesteps << " timesteps in range " << starttime << " to " << endtime << " sec:" << std::endl;
    std::vector<std::vector<vec>> output(2, std::vector<vec>(outputsize));
    output[0][0] = x; output[1][0] = dtx;
    
    // We already have everything for time step 0 so we start at 1:
    int timestepindex = 1;
    for (double t = starttime + timestep; t <= endtime; t = t + timestep)
    {        
        std::cout << timestepindex << "@" << t << "sec" << std::flush;

        mathop::settime(t);
        
        // Nonlinear loop:
        double relchange = 1; int nlit = 0;
        vec xnext = x, dtxnext = dtx;
        while (relchange > tol && (maxnumnlit <= 0 || nlit < maxnumnlit))
        {
            // Solve all formulations that must be solved at the beginning of the nonlinear loop:
            mathop::solve(tosolvebefore);


            // Make all time derivatives available in the universe:
            universe::xdtxdtdtx = {{xnext},{dtxnext},{}};
            
            vec xtolcalc = xnext;
            
            // Reassemble only the non-constant matrices:
            if (isconstant[1] == false || timestepindex == 1)
            {
                myformulation.generatestiffnessmatrix();
                K = myformulation.K(false, false);
            }
            if (isconstant[2] == false || timestepindex == 1)
            {
                myformulation.generatedampingmatrix();
                C = myformulation.C(false, true);
            }
            if (isconstant[0] == false || timestepindex == 1)
            {
                myformulation.generaterhs();
                rhs = myformulation.rhs();
            }
            else
                rhs.updateconstraints();
            
            // Reuse matrices when possible (including the LU decomposition):
            if (isconstant[1] == false || isconstant[2] == false || timestepindex == 1)
            {
                if(order==1)
					leftmat = C + timestep*K;
				if(order==2)
					leftmat = C + 2.0/3.0*timestep*K;
				if(order==3)
					leftmat = C + 6.0/11.0*timestep*K;
				if(order==4)
					leftmat = C + 12.0/25.0*timestep*K;
				if(order==5)
					leftmat = C + 60.0/137.0*timestep*K;
				if(order==6)
					leftmat = C + 60.0/147.0*timestep*K;
                leftmat.reuselu();
            }
            
			//Evaluate the right member
			
			if(order==1)
				right_mbr = C*solution_pool[0] + timestep*rhs;
			if(order==2)
				right_mbr = C*(4.0/3.0*solution_pool[1]-1.0/3.0*solution_pool[0]) + 2.0/3.0*timestep*rhs;
			if(order==3)
				right_mbr = C*(18.0/11.0*solution_pool[2]-9.0/11.0*solution_pool[1]+2.0/11.0*solution_pool[0]) + 6.0/11.0*timestep*rhs;
			if(order==4)
				right_mbr = C*(48.0/25.0*solution_pool[3]-36.0/25.0*solution_pool[2]+16.0/25.0*solution_pool[1]-3.0/25.0*solution_pool[0]) + 12.0/25.0*timestep*rhs;
			if(order==5)
				right_mbr = C*(300.0/137.0*solution_pool[4]-300.0/137.0*solution_pool[3]+200.0/137.0*solution_pool[2]-75.0/137.0*solution_pool[1]+12.0/137.0*solution_pool[0]) + 60.0/137.0*timestep*rhs;
			if(order==6)
				right_mbr = C*(360.0/147.0*solution_pool[5]-450.0/147.0*solution_pool[4]+400.0/147.0*solution_pool[3]-225.0/147.0*solution_pool[2]+72.0/147.0*solution_pool[1]-10.0/147.0*solution_pool[0]) + 60.0/147.0*timestep*rhs;
			
            // Update the solution xnext.
            //xnext = relaxationfactor * mathop::solve(leftmat, right_mbr) + (1.0-relaxationfactor)*solution_pool[order-1];
			xnext = mathop::solve(leftmat, right_mbr) ;
            dtxnext = 1.0/timestep*(xnext-solution_pool[order-1]);            
						
            // Update all fields in the formulation:
            for (int i = 0; i < allfields.size(); i++)
                allfields[i]->setdata(-1, xnext|field(allfields[i]));
            
            relchange = (xnext-xtolcalc).norm()/xnext.norm();
            
            if (islinear == false && verbosity > 0)
                std::cout << " " << relchange << std::flush;

            nlit++; 
            
            
            // Solve all formulations that must be solved at the end of the nonlinear loop:
            mathop::solve(tosolveafter);
            
            
            if (islinear)
                break;
        }
        
		//update the solution pool
		for(int i=0;i<order-1;i++)
			solution_pool[i]=solution_pool[i+1];
		solution_pool[order-1]=xnext;
		
        x = solution_pool[order-1];
        dtx = dtxnext;
        
        if (islinear == false)
            std::cout << " (" << nlit << "NL it)" << std::flush;
        if (timestepindex < numtimesteps-1)
        std::cout << " -> " << std::flush;
        
        // Only one every 'outputeverynthtimestep' solutions is output:
        if (timestepindex%outputeverynthtimestep == 0)
        {
            output[0][timestepindex/outputeverynthtimestep] = x;
            output[1][timestepindex/outputeverynthtimestep] = dtx;
        }
        timestepindex++;
    }
    std::cout << std::endl;
    
    // Remove all time derivatives from the universe:
    universe::xdtxdtdtx = {{},{},{}};
    
    return output;
}


