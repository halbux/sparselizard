#include "jacobian.h"


jacobian::jacobian(elementselector& elemselect, std::vector<double> evaluationcoordinates, expression* meshdeform)
{
    rawfield rf;
    std::vector<densematrix> calced = rf.getjacterms(elemselect, evaluationcoordinates);
    
    field x("x");

    int elementdimension = elemselect.getelementdimension();
    int numberofelements = elemselect.countinselection();
    int problemdimension = universe::mymesh->getmeshdimension();
    int numberofgausspoints = evaluationcoordinates.size()/3;

    double *jac11, *jac12, *jac13, *jac21, *jac22, *jac23, *jac31, *jac32, *jac33, *detjacval;

    
    // Computing JacPoint:
    if (elementdimension == problemdimension - 3)
    {
        switch (problemdimension)
        {
            case 3:

                // This case corresponds to point elements, for which the reference and physical 
                // elements are always the same. The Jacobian determinant must be 1.

                jac[3*0+0] = densematrix(numberofelements,numberofgausspoints,1);
                jac[3*0+1] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*0+2] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*1+0] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*1+1] = densematrix(numberofelements,numberofgausspoints,1);
                jac[3*1+2] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*2+0] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*2+1] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*2+2] = densematrix(numberofelements,numberofgausspoints,1);
                break;
        }
    }
    
    // Computing JacLin:
    if (elementdimension == problemdimension - 2)
    {
        switch (problemdimension)
        {
            case 2:

                // This case corresponds to point elements, for which the reference and physical 
                // elements are always the same. The Jacobian determinant must be 1.

                jac[3*0+0] = densematrix(numberofelements,numberofgausspoints,1);
                jac[3*0+1] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*1+0] = densematrix(numberofelements,numberofgausspoints,0);
                jac[3*1+1] = densematrix(numberofelements,numberofgausspoints,1);
                break;

            case 3:

                // In this case the first row of the jac exists but the second and
                // third ones are equal to zero. In order to be able to compute the determinant
                // and the inverse elements as usual we can compute the second
                // row as a normed vector perpendicular to the first row. The third
                // row will be computed as the cross product of row 1 and row 2 so that row 3 is
                // perpendicular to the first two rows. We then norm row 3 and we are done.
                //
                // Note: Geometrically this corresponds to give to a 1D element
                // a length 1 in dimension 2 and a length 1 in dimension 3. 
                // This does not affect the determinant value, which is just multiplied by 1 x 1.

                jac[3*0+0] = calced[0]; 
                jac[3*0+1] = calced[1];
                jac[3*0+2] = calced[2];
                jac[3*1+0] = densematrix(numberofelements,numberofgausspoints);
                jac[3*1+1] = densematrix(numberofelements,numberofgausspoints);
                jac[3*1+2] = densematrix(numberofelements,numberofgausspoints);
                jac[3*2+0] = densematrix(numberofelements,numberofgausspoints);
                jac[3*2+1] = densematrix(numberofelements,numberofgausspoints);
                jac[3*2+2] = densematrix(numberofelements,numberofgausspoints);
                
                if (meshdeform != NULL)
                {
                    jac[3*0+0].add((meshdeform->getoperationinarray(0,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+1].add((meshdeform->getoperationinarray(1,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+2].add((meshdeform->getoperationinarray(2,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                }

                jac11 = jac[3*0+0].getvalues();
                jac12 = jac[3*0+1].getvalues();
                jac13 = jac[3*0+2].getvalues();
                jac21 = jac[3*1+0].getvalues();
                jac22 = jac[3*1+1].getvalues();
                jac23 = jac[3*1+2].getvalues();
                jac31 = jac[3*2+0].getvalues();
                jac32 = jac[3*2+1].getvalues();
                jac33 = jac[3*2+2].getvalues();
      
                for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                {
                    // Second row perpendicular to first row --> [-b a 0] for row 1 = [a b c]
                    // In case the second row computed with [-b a 0] has a smaller norm than [0 -c b] the latter is used instead.
                    if ( std::sqrt(jac11[i]*jac11[i] + jac12[i]*jac12[i]) > std::sqrt(jac12[i]*jac12[i] + jac13[i]*jac13[i]) )
                    {
                        jac21[i] = -jac12[i];
                        jac22[i] =  jac11[i];
                        jac23[i] =  0;
		    }
                    else
                    {
                        jac21[i] = 0;
                        jac22[i] = -jac13[i];
                        jac23[i] =  jac12[i];
                    }
                    
                    // Norm the second row elements:
                    double normvalue = sqrt(jac21[i]*jac21[i] + jac22[i]*jac22[i] + jac23[i]*jac23[i]);
                    
                    jac21[i] = jac21[i] / normvalue;
                    jac22[i] = jac22[i] / normvalue;
                    jac23[i] = jac23[i] / normvalue;

                    // Third row - normed cross product of row 1 and 2:
                    jac31[i] = jac12[i] * jac23[i] - jac22[i] * jac13[i];
                    jac32[i] = jac13[i] * jac21[i] - jac23[i] * jac11[i];
                    jac33[i] = jac11[i] * jac22[i] - jac21[i] * jac12[i];

                    // Norm the third row elements:
                    normvalue = sqrt(jac31[i]*jac31[i] + jac32[i]*jac32[i] + jac33[i]*jac33[i]);
                    jac31[i] = jac31[i] / normvalue;
                    jac32[i] = jac32[i] / normvalue;
                    jac33[i] = jac33[i] / normvalue;
                }
                break;
        }
    }

    // Computing JacSur:
    if (elementdimension == problemdimension - 1)
    {
        switch (problemdimension)
        {
            case 1:

                // This case corresponds to point elements, for which the reference and physical 
                // elements are always the same. The Jacobian determinant must be 1.
                jac[3*0+0] = densematrix(numberofelements,numberofgausspoints,1);
                break;

            case 2:
                
                // In this case the first row of the jac exists but the second one
                // is equal to zero. In order to be able to compute the determinant
                // and the inverse elements as usual we can compute the second
                // row as the normed vector perpendicular to the first row (for that we
                // know that [-b a] is perpendicular to any vector [a b]).
                //
                // Note: Geometrically this corresponds to give to a 
                // 1D element a length 1 in dimension 2. 
                // This does not affect the determinant value, which is just multiplied by 1.
                    
                jac[3*0+0] = calced[0]; 
                jac[3*0+1] = calced[1];
                jac[3*1+0] = densematrix(numberofelements,numberofgausspoints);
                jac[3*1+1] = densematrix(numberofelements,numberofgausspoints);

                if (meshdeform != NULL)
                {
                    jac[3*0+0].add((meshdeform->getoperationinarray(0,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+1].add((meshdeform->getoperationinarray(1,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                }
                
                jac11 = jac[3*0+0].getvalues();
                jac12 = jac[3*0+1].getvalues();
                jac21 = jac[3*1+0].getvalues();
                jac22 = jac[3*1+1].getvalues();

                for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                {
                    // Second row perpendicular to first row --> [-b a] for row 1 = [a b]:
                    jac21[i] = -jac12[i];
                    jac22[i] = jac11[i];

                    // Norm the second row elements:
                    double normvalue = sqrt(jac21[i]*jac21[i] + jac22[i]*jac22[i]);
                    jac21[i] = jac21[i] / normvalue;
                    jac22[i] = jac22[i] / normvalue;    
                }
                break;

            case 3:
    
                // In this case the first 2 rows of the jac exist but the third
                // one is equal to zero. In order to be able to compute the determinant
                // and the inverse elements as usual we can compute the third
                // row as the cross product of row 1 and row 2 so that row 3 is
                // perpendicular to the first two rows. We then norm row 3 and we are done.
                //
                // Note: Geometrically this corresponds to give to a 2D element
                // a thickness of 1 in the third dimension which does not affect the
                // determinant value which is just multiplied by 1.
                
                jac[3*0+0] = calced[0]; 
                jac[3*0+1] = calced[1];
                jac[3*0+2] = calced[2];
                jac[3*1+0] = calced[3];
                jac[3*1+1] = calced[4];
                jac[3*1+2] = calced[5];
                jac[3*2+0] = densematrix(numberofelements,numberofgausspoints);
                jac[3*2+1] = densematrix(numberofelements,numberofgausspoints);
                jac[3*2+2] = densematrix(numberofelements,numberofgausspoints);
                
                if (meshdeform != NULL)
                {
                    jac[3*0+0].add((meshdeform->getoperationinarray(0,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+1].add((meshdeform->getoperationinarray(1,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+2].add((meshdeform->getoperationinarray(2,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+0].add((meshdeform->getoperationinarray(0,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+1].add((meshdeform->getoperationinarray(1,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+2].add((meshdeform->getoperationinarray(2,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                }
                
                jac11 = jac[3*0+0].getvalues();
                jac12 = jac[3*0+1].getvalues();
                jac13 = jac[3*0+2].getvalues();
                jac21 = jac[3*1+0].getvalues();
                jac22 = jac[3*1+1].getvalues();
                jac23 = jac[3*1+2].getvalues();
                jac31 = jac[3*2+0].getvalues();
                jac32 = jac[3*2+1].getvalues();
                jac33 = jac[3*2+2].getvalues();
                
                for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                {
                    // Third row - normed cross product of row 1 and 2:
                    jac31[i] = jac12[i] * jac23[i] - jac22[i] * jac13[i];
                    jac32[i] = jac13[i] * jac21[i] - jac23[i] * jac11[i];
                    jac33[i] = jac11[i] * jac22[i] - jac21[i] * jac12[i];

                    // Norm the third row elements:
                    double normvalue = sqrt(jac31[i]*jac31[i] + jac32[i]*jac32[i] + jac33[i]*jac33[i]);
                    jac31[i] = jac31[i] / normvalue;
                    jac32[i] = jac32[i] / normvalue;
                    jac33[i] = jac33[i] / normvalue;      
                }
                break;
        }
    }

    // Computing JacVol:
    if (elementdimension == problemdimension)
    {
        switch (problemdimension)
        {
            case 1:

                jac[3*0+0] = calced[0]; 
                
                if (meshdeform != NULL)
                    jac[3*0+0].add((meshdeform->getoperationinarray(0,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                break; 
                
            case 2:

                jac[3*0+0] = calced[0]; 
                jac[3*0+1] = calced[1];
                jac[3*1+0] = calced[2];
                jac[3*1+1] = calced[3];
                
                if (meshdeform != NULL)
                {
                    jac[3*0+0].add((meshdeform->getoperationinarray(0,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+1].add((meshdeform->getoperationinarray(1,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+0].add((meshdeform->getoperationinarray(0,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+1].add((meshdeform->getoperationinarray(1,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                }
                break;
                
            case 3:

                jac[3*0+0] = calced[0]; 
                jac[3*0+1] = calced[1];
                jac[3*0+2] = calced[2];
                jac[3*1+0] = calced[3];
                jac[3*1+1] = calced[4]; 
                jac[3*1+2] = calced[5];
                jac[3*2+0] = calced[6];
                jac[3*2+1] = calced[7];
                jac[3*2+2] = calced[8]; 
                
                if (meshdeform != NULL)
                {
                    jac[3*0+0].add((meshdeform->getoperationinarray(0,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+1].add((meshdeform->getoperationinarray(1,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*0+2].add((meshdeform->getoperationinarray(2,0)->interpolate(1, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+0].add((meshdeform->getoperationinarray(0,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+1].add((meshdeform->getoperationinarray(1,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*1+2].add((meshdeform->getoperationinarray(2,0)->interpolate(2, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*2+0].add((meshdeform->getoperationinarray(0,0)->interpolate(3, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*2+1].add((meshdeform->getoperationinarray(1,0)->interpolate(3, elemselect, evaluationcoordinates))[1][0]);
                    jac[3*2+2].add((meshdeform->getoperationinarray(2,0)->interpolate(3, elemselect, evaluationcoordinates))[1][0]);
                }
                break;
        }
    }
    
    
    
    // Compute the determinant:
    detjac = densematrix(numberofelements,numberofgausspoints);

    switch (problemdimension)
    {
        case 1:            
            
            detjac = jac[3*0+0].copy();
            break;
            
        case 2:

            jac11 = jac[3*0+0].getvalues();
            jac12 = jac[3*0+1].getvalues();
            jac21 = jac[3*1+0].getvalues();
            jac22 = jac[3*1+1].getvalues();
            
            detjacval = detjac.getvalues();
        
            for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                detjacval[i] = jac11[i] * jac22[i] - jac12[i] * jac21[i];
            break;
            
        case 3:

            jac11 = jac[3*0+0].getvalues();
            jac12 = jac[3*0+1].getvalues();
            jac13 = jac[3*0+2].getvalues();
            jac21 = jac[3*1+0].getvalues();
            jac22 = jac[3*1+1].getvalues();
            jac23 = jac[3*1+2].getvalues();
            jac31 = jac[3*2+0].getvalues();
            jac32 = jac[3*2+1].getvalues();
            jac33 = jac[3*2+2].getvalues();
            
            detjacval = detjac.getvalues();

            for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                detjacval[i] = jac11[i] * (jac22[i] * jac33[i] - jac32[i] * jac23[i]) - jac12[i] * (jac21[i] * jac33[i] - jac31[i] * jac23[i]) + jac13[i] * (jac21[i] * jac32[i] - jac31[i] * jac22[i]);
            break;
            
    }

    if (universe::isaxisymmetric)
    {
        xcoord = (x.getpointer()->interpolate(0, 0, elemselect, evaluationcoordinates))[1][0];
        if (meshdeform != NULL)
            xcoord.add((meshdeform->getoperationinarray(0,0)->interpolate(0, elemselect, evaluationcoordinates))[1][0]);
        jac[3*2+2] = xcoord.copy();
    }
}

jacobian jacobian::extractsubset(std::vector<int>& selectedelementindexes)
{
    jacobian subjac;

    subjac.detjac = detjac.extractrows(selectedelementindexes);
    if (xcoord.isdefined())
        subjac.xcoord = xcoord.extractrows(selectedelementindexes);

    for (int i = 0; i < jac.size(); i++)
    {
        if (jac[i].isdefined())
            subjac.jac[i] = jac[i].extractrows(selectedelementindexes);
    }

    subjac.invjac = std::vector<densematrix>(invjac.size());
    for (int i = 0; i < invjac.size(); i++)
    {
        if (invjac[i].isdefined())
            subjac.invjac[i] = invjac[i].extractrows(selectedelementindexes);
    }

    return subjac;
}

densematrix jacobian::getdetjac(void)
{ 
    densematrix detj = detjac.copy();

    if (universe::isaxisymmetric)
        detj.multiplyelementwise(xcoord);

    return detj;
}    

densematrix jacobian::getjac(int row, int column) { return jac[3*row+column]; }

densematrix jacobian::getinvjac(int row, int column)
{
    // If the Jacobian inverse has not been computed yet compute it:
    if (invjac.size() == 0)
    {
        invjac = std::vector<densematrix>(3*3);

        int numberofelements = detjac.countrows();
        int numberofgausspoints = detjac.countcolumns();

        double* jac11, *jac12, *jac13, *jac21, *jac22, *jac23, *jac31, *jac32, *jac33, *detjacval;
        double* invjac11, *invjac12, *invjac13, *invjac21, *invjac22, *invjac23, *invjac31, *invjac32, *invjac33;

        int problemdimension = universe::mymesh->getmeshdimension();
        
        switch (problemdimension)
        {
            case 1:
                invjac[3*0+0] = densematrix(numberofelements,numberofgausspoints);

                jac11 = jac[3*0+0].getvalues();
                invjac11 = invjac[3*0+0].getvalues();

                for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                    invjac11[i] = 1.0/jac11[i];
                break;
            case 2:
                invjac[3*0+0] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*0+1] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*1+0] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*1+1] = densematrix(numberofelements,numberofgausspoints);

                jac11 = jac[3*0+0].getvalues();
                jac12 = jac[3*0+1].getvalues();
                jac21 = jac[3*1+0].getvalues();
                jac22 = jac[3*1+1].getvalues();

                invjac11 = invjac[3*0+0].getvalues();
                invjac12 = invjac[3*0+1].getvalues();
                invjac21 = invjac[3*1+0].getvalues();
                invjac22 = invjac[3*1+1].getvalues();

                detjacval = detjac.getvalues();

                for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                {
                    invjac11[i] =   jac22[i] / detjacval[i];
                    invjac12[i] = - jac12[i] / detjacval[i];
                    invjac21[i] = - jac21[i] / detjacval[i];
                    invjac22[i] =   jac11[i] / detjacval[i];
                }
                break;
            case 3:
                invjac[3*0+0] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*0+1] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*0+2] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*1+0] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*1+1] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*1+2] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*2+0] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*2+1] = densematrix(numberofelements,numberofgausspoints);
                invjac[3*2+2] = densematrix(numberofelements,numberofgausspoints);

                jac11 = jac[3*0+0].getvalues();
                jac12 = jac[3*0+1].getvalues();
                jac13 = jac[3*0+2].getvalues();
                jac21 = jac[3*1+0].getvalues();
                jac22 = jac[3*1+1].getvalues();
                jac23 = jac[3*1+2].getvalues();
                jac31 = jac[3*2+0].getvalues();
                jac32 = jac[3*2+1].getvalues();
                jac33 = jac[3*2+2].getvalues();

                invjac11 = invjac[3*0+0].getvalues();
                invjac12 = invjac[3*0+1].getvalues();
                invjac13 = invjac[3*0+2].getvalues();
                invjac21 = invjac[3*1+0].getvalues();
                invjac22 = invjac[3*1+1].getvalues();
                invjac23 = invjac[3*1+2].getvalues();
                invjac31 = invjac[3*2+0].getvalues();
                invjac32 = invjac[3*2+1].getvalues();
                invjac33 = invjac[3*2+2].getvalues();

                detjacval = detjac.getvalues();

                for (int i = 0; i < numberofelements*numberofgausspoints; i++)
                {
                    invjac11[i] = (jac22[i] * jac33[i] - jac23[i] * jac32[i]) /detjacval[i];
                    invjac12[i] = (jac13[i] * jac32[i] - jac12[i] * jac33[i]) /detjacval[i];
                    invjac13[i] = (jac12[i] * jac23[i] - jac13[i] * jac22[i]) /detjacval[i];
                    invjac21[i] = (jac23[i] * jac31[i] - jac21[i] * jac33[i]) /detjacval[i];
                    invjac22[i] = (jac11[i] * jac33[i] - jac13[i] * jac31[i]) /detjacval[i];
                    invjac23[i] = (jac13[i] * jac21[i] - jac11[i] * jac23[i]) /detjacval[i];
                    invjac31[i] = (jac21[i] * jac32[i] - jac22[i] * jac31[i]) /detjacval[i];
                    invjac32[i] = (jac12[i] * jac31[i] - jac11[i] * jac32[i]) /detjacval[i];
                    invjac33[i] = (jac11[i] * jac22[i] - jac12[i] * jac21[i]) /detjacval[i];
                }
                break;
        }
        
        if (universe::isaxisymmetric)
        {
            invjac[3*2+2] = jac[3*2+2].copy();
            invjac[3*2+2].invert();
        }
    }

    return invjac[3*row+column];
}





