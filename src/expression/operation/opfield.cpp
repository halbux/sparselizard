#include "opfield.h"


void opfield::setspacederivative(int whichderivative)
{
    // Make sure a single space derivative is applied.
    if (spacederivative != 0 || kietaphiderivative != 0)
    {
        std::cout << "Error in 'opfield' object: cannot apply more than one space derivative to a field" << std::endl;
        abort();
    }
    spacederivative = whichderivative;
}

void opfield::setkietaphiderivative(int whichderivative)
{
    // Make sure a single space derivative is applied.
    if (spacederivative != 0 || kietaphiderivative != 0)
    {
        std::cout << "Error in 'opfield' object: cannot apply more than one space derivative to a field" << std::endl;
        abort();
    }
    kietaphiderivative = whichderivative;
}

void opfield::increasetimederivativeorder(int amount)
{
    timederivativeorder += amount;

    if (not(myfield->ismultiharmonic()) && timederivativeorder > 2)
    {
        std::cout << "Error in 'opfield' object: time derivative order can exceed 2 only for multiharmonic fields" << std::endl;
        abort();
    }
}

bool opfield::isharmonicone(std::vector<int> disjregs)
{
    std::vector<int> myharms = myfield->getharmonics();
    return (myharms.size() == 1 && myharms[0] == 1);
}

std::vector<std::vector<densematrix>> opfield::interpolate(int kietaphiderivative, elementselector& elemselect, std::vector<double>& evaluationcoordinates)
{
    if (timederivativeorder == 0 && spacederivative == 0)
        return myfield->interpolate(kietaphiderivative, formfunctioncomponent, elemselect, evaluationcoordinates);
    else
    {
        std::cout << "Error in 'opfield' object: expression provided for mesh deformation is invalid" << std::endl;
        std::cout << "Operation was:" << std::endl;
        this->print();
        abort();
    }
}

std::vector<std::vector<densematrix>> opfield::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }

    // Time derivative of non-multiharmonic fields:
    std::shared_ptr<coefmanager> cmbkp = NULL;
    if (myfield->ismultiharmonic() == false && timederivativeorder > 0)
    {
        // Get the vector in the universe corresponding to the field time derivative.
        // This was created and set to the universe by a time resolution object.
        if ((universe::xdtxdtdtx)[timederivativeorder].size() == 0)
        {
            std::vector<std::string> messtr = {"","dt","dtdt"};
            std::cout << "Error in 'opfield' object: the " << messtr[timederivativeorder] << "(";
            myfield->print();
            std::cout << ") value was not made available by a time resolution or by a call to 'settimederivative'" << std::endl;
            abort();
        }
        cmbkp = myfield->harmonic(1)->resetcoefmanager();
        // Set the field value to the field time derivative value on all regions:
        myfield->setdata(-1, (universe::xdtxdtdtx)[timederivativeorder][0]|field(myfield));
    }

    std::vector<std::vector<densematrix>> output;

    // In case there is no space derivative applied:
    if (spacederivative == 0 && kietaphiderivative == 0)
        output = myfield->interpolate(0, formfunctioncomponent, elemselect, evaluationcoordinates);
    else
    {
        if (kietaphiderivative != 0)
            output = myfield->interpolate(kietaphiderivative, formfunctioncomponent, elemselect, evaluationcoordinates);
        else
        {
            // Otherwise compute the x, y or z field derivative using ki, eta,
            // phi derivatives in the reference element and invjac() terms.

            // Compute the Jacobian terms or reuse if available in the universe.
            std::shared_ptr<jacobian> myjac;
            if (universe::isreuseallowed && universe::computedjacobian != NULL)
                myjac = universe::computedjacobian;
            else
                myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));

            if (universe::isreuseallowed)
                universe::computedjacobian = myjac;

            // Compute the required ki, eta and phi derivatives:
            std::vector<std::vector<densematrix>> dkiargmat, detaargmat, dphiargmat;

            // Get the element dimension in the selected elements:
            int elementdimension = elemselect.getelementdimension();

                dkiargmat = myfield->interpolate(1, formfunctioncomponent, elemselect, evaluationcoordinates);
            if (elementdimension > 1)
                detaargmat = myfield->interpolate(2, formfunctioncomponent, elemselect, evaluationcoordinates);
            if (elementdimension > 2)
                dphiargmat = myfield->interpolate(3, formfunctioncomponent, elemselect, evaluationcoordinates);

            // Computing the x, y and z derivatives for all harmonics based on the
            // Jacobian matrix and the computed ki, eta and phi derivatives.
            // Skip the sin0 harmonic.
            for (int harm = 1; harm < dkiargmat.size(); harm++)
            {
                // Skip any non existent harmonic:
                if (dkiargmat[harm].size() == 0)
                    continue;

                    dkiargmat[harm][0].multiplyelementwise(myjac->getinvjac(spacederivative-1,0));
                if (elementdimension > 1)
                {
                    detaargmat[harm][0].multiplyelementwise(myjac->getinvjac(spacederivative-1,1));
                    dkiargmat[harm][0].add(detaargmat[harm][0]);
                }
                if (elementdimension > 2)
                {
                    dphiargmat[harm][0].multiplyelementwise(myjac->getinvjac(spacederivative-1,2));
                    dkiargmat[harm][0].add(dphiargmat[harm][0]);
                }
            }
            output = dkiargmat;
        }
    }

    if (cmbkp != NULL)
        myfield->harmonic(1)->setcoefmanager(cmbkp);

    if (myfield->ismultiharmonic() && timederivativeorder > 0)
        output = harmonic::timederivative(timederivativeorder, output);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), output);

    return output;
}

densematrix opfield::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }

    std::vector<std::vector<densematrix>> interpolatedfield = interpolate(elemselect, evaluationcoordinates, meshdeform);
    // Compute at 'numtimevals' instants in time the multiharmonic field:
    densematrix output = myfft::inversefft(interpolatedfield, numtimeevals, elemselect.countinselection(), evaluationcoordinates.size()/3);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    return output;
}

bool opfield::isvalueorientationdependent(std::vector<int> disjregs)
{
    if (myfield->gettypename() == "x" || myfield->gettypename() == "y" || myfield->gettypename() == "z")
        return false;

    for (int i = 0; i < disjregs.size(); i++)
    {
        int elementtypenumber = (universe::mymesh->getdisjointregions())->getelementtypenumber(disjregs[i]);
        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, myfield->gettypename());
        if ( myformfunction->isorientationdependent(myfield->getinterpolationorder(disjregs[i])) )
            return true;
    }
    return false;
}

std::shared_ptr<operation> opfield::copy(void)
{
    std::shared_ptr<opfield> op(new opfield(myfield));
    *op = *this;
    return op;
}

std::vector<double> opfield::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::string mytype = myfield->gettypename();
    if (timederivativeorder != 0 || spacederivative != 0 || kietaphiderivative != 0)
    {
        std::cout << "Error in 'opfield' object: evaluate does not allow derivatives" << std::endl;
        abort();
    }
    if (mytype == "x")
        return xcoords;
    if (mytype == "y")
        return ycoords;
    if (mytype == "z")
        return zcoords;

    std::cout << "Error in 'opfield' object: evaluate only allows the x, y and z field" << std::endl;
    abort();
}

void opfield::print(void)
{
    for (int i = 0; i < timederivativeorder; i++)
        std::cout << "dt";

    std::vector<std::string> nonedxdydz = {"","dx","dy","dz"};
    std::vector<std::string> nonecompxcompycompz = {"compx","compy","compz"};

    std::cout << nonedxdydz[spacederivative];
    // For fields without subfields:
    if (fieldcomponent == -1)
    {
        if (myfield->countformfunctioncomponents() > 1)
            std::cout << nonecompxcompycompz[formfunctioncomponent];
    }
    else
        std::cout << nonecompxcompycompz[fieldcomponent];

    myfield->print();
}
