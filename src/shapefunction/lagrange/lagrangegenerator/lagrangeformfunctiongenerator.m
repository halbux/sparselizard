function lagrangeformfunctiongenerator(elementname, elementtypenumber, order, coordinateorder)

% This function writes all Lagrange form functions for element 'elementname'
% up to order 'order' in 'lagrangeELEMENTNAME.cpp' and '.h'.

elementname = lower(elementname);
filename = ['lagrange' elementname];


%%%%% Write first the header:

fileid = fopen([filename '.h'], 'W');

fprintf(fileid, '#ifndef %s_H\n', upper(filename));
fprintf(fileid, '#define %s_H\n', upper(filename));
fprintf(fileid, '\n');
fprintf(fileid, '#include <vector>\n');
fprintf(fileid, '#include "polynomial.h"\n');
fprintf(fileid, '\n');
fprintf(fileid, '\n');
fprintf(fileid, 'namespace %s\n', filename);
fprintf(fileid, '{\n');
% fprintf(fileid, '    // ''count(order)'' gives the number of Lagrange form functions at ''order'':\n');
% fprintf(fileid, '    int count(int);\n');
fprintf(fileid, '    // ''getnodecoordinates'' gives the ki, eta and phi \n    // coordinates of every Lagrange node for the required order:\n');
fprintf(fileid, '    std::vector<double> getnodecoordinates(int order);\n');
fprintf(fileid, '    // ''getformfunctionpolynomials'' outputs a vector giving the form function polynomials:\n');
fprintf(fileid, '    std::vector<polynomial> getformfunctionpolynomials(int order);\n');
fprintf(fileid, '}\n');
fprintf(fileid, '\n');
fprintf(fileid, '#endif');

fclose(fileid);


%%%%% Write the first part of the .cpp file:

fileid = fopen([filename '.cpp'], 'W');

fprintf(fileid, '#include <vector>\n');
fprintf(fileid, '#include "math.h"\n');
fprintf(fileid, '#include "element.h"\n');
fprintf(fileid, '#include "polynomial.h"\n');
fprintf(fileid, '#include "%s.h"\n', filename);
fprintf(fileid, '\n');
fprintf(fileid, '\n');


% %%%%% Write the 'count' function:
% 
% fprintf(fileid, 'int %s::count(int order)\n', filename);
% fprintf(fileid, '{\n');
% % Compute the number of Lagrange form functions for an element type and order:
% if (strcmp(elementname,'point'))
%     numlagrangeformfunc = '1';
% end
% if (strcmp(elementname,'line'))
%     numlagrangeformfunc = 'order+1';
% end
% if (strcmp(elementname,'triangle'))
%     numlagrangeformfunc = '0.5*( pow(order+1,2) + (order+1) )';
% end
% if (strcmp(elementname,'quadrangle'))
%     numlagrangeformfunc = 'pow(order+1,2)';
% end
% if (strcmp(elementname,'tetrahedron'))
%     numlagrangeformfunc = '1.0/6.0*(order+1)*(order+2)*(order+3)';
% end
% if (strcmp(elementname,'hexahedron'))
%     numlagrangeformfunc = 'pow(order+1,3)';
% end
% if (strcmp(elementname,'prism'))
%     numlagrangeformfunc = '0.5*( pow(order+1,2) + (order+1) ) * (order+1)';
% end
% 
% fprintf(fileid, '     return %s;\n', numlagrangeformfunc);
% fprintf(fileid, '}\n');
% fprintf(fileid, '\n');


%%%%% Write the 'getnodecoordinates' function:


fprintf(fileid, 'std::vector<double> %s::getnodecoordinates(int order)\n', filename);
fprintf(fileid, '{\n');

% Loop on all orders:
fprintf(fileid, '     switch (order)\n');
fprintf(fileid, '     {\n');

for currentorder = 1:coordinateorder
    ['ORDER ' num2str(currentorder)]
    
    fprintf(fileid, '          case %d:\n', currentorder);

    % Default values:
    ki = []; eta = []; phi = [];
    
    if strcmp(elementname,'line')
        
        numberofnodes = currentorder+1;
        
        elementdimension = 1;
        % Corner nodes:
        ki = [-1 1];
        
        ki = getnodesinteriortoedge([ki; zeros(1,2); zeros(1,2)], currentorder, 1);
        
    end
    
    if strcmp(elementname,'triangle')
        
        numberofnodes = 0.5*( (currentorder+1)^2 + (currentorder+1) );
        
        elementdimension = 2;
        % Corner nodes:
        ki =  vpa([0 1 0]);
        eta = vpa([0 0 1]);

        [ki, eta] = getnodesinteriortotriangularface([ki; eta; zeros(1,3)], currentorder, 1);            

    end
    
    if strcmp(elementname,'quadrangle')
        
        numberofnodes = (currentorder+1)^2;
        
        elementdimension = 2;
        % Corner nodes:
        ki =  vpa([-1 1 1 -1]);
        eta = vpa([-1 -1 1 1]);

        [ki, eta] = getnodesinteriortoquadrangularface([ki; eta; zeros(1,4)], currentorder, 1);            

    end
    
    if strcmp(elementname,'tetrahedron')
        
        numberofnodes = 1/6*(currentorder+1)*(currentorder+2)*(currentorder+3);
        
        elementdimension = 3;
        % Corner nodes:
        ki =  vpa([0 1 0 0]);
        eta = vpa([0 0 1 0]);
        phi = vpa([0 0 0 1]);
        
        [ki, eta, phi] = getnodesinteriortotetrahedronvolume([ki; eta; phi], currentorder, 1);            

    end
    
    if strcmp(elementname,'hexahedron')
        
        numberofnodes = (currentorder+1)^3;
        
        elementdimension = 3;
        % Corner nodes:
        ki =  vpa([-1 1 1 -1 -1 1 1 -1]);
        eta = vpa([-1 -1 1 1 -1 -1 1 1]);
        phi = vpa([-1 -1 -1 -1 1 1 1 1]);
        
        [ki, eta, phi] = getnodesinteriortohexahedronvolume([ki; eta; phi], currentorder, 1);            

    end
    
    if strcmp(elementname,'prism')
        
        numberofnodes = 0.5*( (currentorder+1)^2 + (currentorder+1)  *  (currentorder+1) );

        elementdimension = 3;
        % Corner nodes:
        ki =  vpa([0 1 0 0 1 0]);
        eta = vpa([0 0 1 0 0 1]);
        phi = vpa([-1 -1 -1 1 1 1]);
        
        [ki, eta, phi] = getnodesinteriortoprismvolume([ki; eta; phi], currentorder, 1);            

    end
    
    if isempty(eta)
        eta = zeros(1,length(ki));
    end
    if isempty(phi)
        phi = vpa(zeros(1,length(ki)));
    end
    
    fprintf(fileid, '               return std::vector<double> {%s, %s, %s', char(vpa(ki(1),20)), char(vpa(eta(1),20)), char(vpa(phi(1),20)));
    for i = 2:length(ki)
        fprintf(fileid, ', %s, %s, %s', char(vpa(ki(i),20)), char(vpa(eta(i),20)), char(vpa(phi(i),20)));
    end
    fprintf(fileid, '};\n');
    
    'Node coordinates ki eta phi:'
    ki
    eta
    phi
    
end

fprintf(fileid, '          default:\n');

fprintf(fileid, '               std::cout << "Error in ''%s'' namespace: coordinates of order %d and above not defined" << std::endl;\n', filename, coordinateorder+1);
fprintf(fileid, '               abort();\n');
fprintf(fileid, '               break;\n');


fprintf(fileid, '     }\n');
fprintf(fileid, '}\n\n');


fprintf(fileid, 'std::vector<polynomial> %s::getformfunctionpolynomials(int order)\n', filename);
fprintf(fileid, '{\n');


%%%%% Initialise the output vector size:

fprintf(fileid, '     element %s(%d,order);\n', elementname, elementtypenumber);
fprintf(fileid, '     std::vector<polynomial> formfunctionpoly(%s.countcurvednodes());\n', elementname);

fprintf(fileid, '\n');

%%%%% Loop on all orders:

fprintf(fileid, '     switch (order)\n');
fprintf(fileid, '     {\n');

for currentorder = 1:order
    ['ORDER ' num2str(currentorder)]

    fprintf(fileid, '          case %d:\n', currentorder);
    
    %%%%% Get the node coordinates for the Lagrange form functions at the current order:
    
    % Default values:
    ki = []; eta = []; phi = [];
    
    if strcmp(elementname,'line')
        
        numberofnodes = currentorder+1;
        
        elementdimension = 1;
        % Corner nodes:
        ki = [-1 1];
        
        ki = getnodesinteriortoedge([ki; zeros(1,2); zeros(1,2)], currentorder, 1);
        
    end
    
    if strcmp(elementname,'triangle')
        
        numberofnodes = 0.5*( (currentorder+1)^2 + (currentorder+1) );
        
        elementdimension = 2;
        % Corner nodes:
        ki =  vpa([0 1 0]);
        eta = vpa([0 0 1]);

        [ki, eta] = getnodesinteriortotriangularface([ki; eta; zeros(1,3)], currentorder, 1);            

    end
    
    if strcmp(elementname,'quadrangle')
        
        numberofnodes = (currentorder+1)^2;
        
        elementdimension = 2;
        % Corner nodes:
        ki =  vpa([-1 1 1 -1]);
        eta = vpa([-1 -1 1 1]);

        [ki, eta] = getnodesinteriortoquadrangularface([ki; eta; zeros(1,4)], currentorder, 1);            

    end
    
    if strcmp(elementname,'tetrahedron')
        
        numberofnodes = 1/6*(currentorder+1)*(currentorder+2)*(currentorder+3);
        
        elementdimension = 3;
        % Corner nodes:
        ki =  vpa([0 1 0 0]);
        eta = vpa([0 0 1 0]);
        phi = vpa([0 0 0 1]);
        
        [ki, eta, phi] = getnodesinteriortotetrahedronvolume([ki; eta; phi], currentorder, 1);            

    end
    
    if strcmp(elementname,'hexahedron')
        
        numberofnodes = (currentorder+1)^3;
        
        elementdimension = 3;
        % Corner nodes:
        ki =  vpa([-1 1 1 -1 -1 1 1 -1]);
        eta = vpa([-1 -1 1 1 -1 -1 1 1]);
        phi = vpa([-1 -1 -1 -1 1 1 1 1]);
        
        [ki, eta, phi] = getnodesinteriortohexahedronvolume([ki; eta; phi], currentorder, 1);            

    end
    
    if strcmp(elementname,'prism')
        
        numberofnodes = 0.5*( (currentorder+1)^2 + (currentorder+1) ) * (currentorder+1);

        elementdimension = 3;
        % Corner nodes:
        ki =  vpa([0 1 0 0 1 0]);
        eta = vpa([0 0 1 0 0 1]);
        phi = vpa([-1 -1 -1 1 1 1]);
        
        [ki, eta, phi] = getnodesinteriortoprismvolume([ki; eta; phi], currentorder, 1);            

    end
    
    if isempty(eta)
        eta = zeros(1,length(ki));
    end
    if isempty(phi)
        phi = vpa(zeros(1,length(ki)));
    end
    
%     figure;
%     a=0;
%     if elementdimension == 3
%         scatter3(ki,eta,phi, 'filled');
%         nodenum = [1:length(ki)]'-a; nodenumstr = num2str(nodenum); nodelabels = cellstr(nodenumstr);
%         % Displacement so the text does not overlay the data points:
%         dki = 0.02; deta = 0.02; dphi = 0.02;
%         text(double(ki+dki), double(eta+deta), double(phi+dphi), nodelabels);
%     else
%         scatter(ki,eta, 'filled');
%         nodenum = [1:length(ki)]'-a; nodenumstr = num2str(nodenum); nodelabels = cellstr(nodenumstr);
%         % Displacement so the text does not overlay the data points:
%         dki = 0.02; deta = 0.02;
%         text(double(ki+dki), double(eta+deta), nodelabels);
%     end
%     drawnow;


    % The general polynomial is defined by the number of nodes it has to go go through

    formfunctionmonomialvariablematrix = [];
    for i = (0:currentorder)
        formfunctionmonomialvariablematrix = [formfunctionmonomialvariablematrix; getallpossiblevectorswithsumequaln(elementdimension, i)];
    end
    
    if strcmp(elementname, 'quadrangle')
        % Must go till 2*currentorder:
        for i = (currentorder + 1:2*currentorder)
            formfunctionmonomialvariablematrix = [formfunctionmonomialvariablematrix; getallpossiblevectorswithsumequaln(elementdimension, i)];
        end
        % Remove the full > currentorder1 monomials:
        formfunctionmonomialvariablematrix(any(formfunctionmonomialvariablematrix > currentorder, 2),:) = [];
    end
    
    if elementdimension == 1
        formfunctionmonomialvariablematrix(:,2:3) = 0;
    end
    if elementdimension == 2
        formfunctionmonomialvariablematrix(:,3) = 0;
    end
    
    if strcmp(elementname, 'hexahedron')
        % Must go till 3*currentorder:
        for i = (currentorder + 1:3*currentorder)
            formfunctionmonomialvariablematrix = [formfunctionmonomialvariablematrix; getallpossiblevectorswithsumequaln(elementdimension, i)];
        end
        % Remove the full > currentorder1 monomials:
        formfunctionmonomialvariablematrix(any(formfunctionmonomialvariablematrix > currentorder, 2),:) = [];
    end
    if strcmp(elementname, 'prism')
        % Must take the triangle form functions (linear without phi) then x phi, phi^2, ... to the power required by the z direction edges:
        % First keep only the terms with 0 order for phi:
        formfunctionmonomialvariablematrix(formfunctionmonomialvariablematrix(:,3) > 0,:) = [];
        % Now add the products by phi, phi^2,...:
        formfunctionmonomialvariablematrixphiorder0 = formfunctionmonomialvariablematrix;
        for phiorder = 1:currentorder
            formfunctionmonomialvariablematrix = [formfunctionmonomialvariablematrix; formfunctionmonomialvariablematrixphiorder0(:,1:2) formfunctionmonomialvariablematrixphiorder0(:,3)+phiorder];
        end
    end
    
    % 'monomialvariables' must have a row per node:
    if size(formfunctionmonomialvariablematrix, 1) ~= numberofnodes
        error('There are not as many monomials as nodes');
    end

    % Now 'monomialvariables' column 1 gives the 'ki' exponent for every
    % monomial, column is for eta and 3 for phi, if needed


    %%%%%%%%%% Creating the Vandermonde matrix
    % Construction is explained here: http://femwiki.wikidot.com/elements:lagrange-elements

    vandermondematrix = vpa(zeros(numberofnodes, numberofnodes));

    % Evaluate every monomial (Vandermonde column) at every node (Vandermonde row):
    for i = (1:numberofnodes)

        if elementdimension == 1
            vandermondematrix(:,i) = transpose( ki.^formfunctionmonomialvariablematrix(i,1) );
        end
        if elementdimension == 2
            vandermondematrix(:,i) = transpose( ki.^formfunctionmonomialvariablematrix(i,1).*eta.^formfunctionmonomialvariablematrix(i,2) );
        end
        if elementdimension == 3
            vandermondematrix(:,i) = transpose( ki.^formfunctionmonomialvariablematrix(i,1).*eta.^formfunctionmonomialvariablematrix(i,2).*phi.^formfunctionmonomialvariablematrix(i,3) );
        end

    end


    %%%%%%%%%% Get the polynomial coefficients for every lagrange form function
    % The coefficients for e.g. form function 2 are Vandermonde\[0;1;0;0;0...0]

    formfunctioncoefficientmatrix = vpa(zeros(numberofnodes, numberofnodes));

    for i = (1:numberofnodes)

        rhs = zeros(numberofnodes,1);
        rhs(i) = 1;
        formfunctioncoefficientmatrix(i,:) = transpose(vandermondematrix\rhs);

    end
    
    formfunctioncoefficientmatrix(abs(double(formfunctioncoefficientmatrix)) < 1e-20) = 0;

    formfunctioncoefficientmatrix
    

    maxkiorder = max(formfunctionmonomialvariablematrix(:,1));
    maxetaorder = max(formfunctionmonomialvariablematrix(:,2));
    maxphiorder = max(formfunctionmonomialvariablematrix(:,3));

    formfunctionmonomialvariablematrix
    for i = 1:numberofnodes
        fprintf(fileid, '               formfunctionpoly[%d].set({', i-1);
    
        for currentkiorder = 0:maxkiorder
            
            if isempty(find(formfunctionmonomialvariablematrix(:,1) == currentkiorder, 1))
                continue;
            end
            
            if currentkiorder > 0
                fprintf(fileid, ', ');
            end
            fprintf(fileid, '{');
            for currentetaorder = 0:maxetaorder
                
                if isempty(find(formfunctionmonomialvariablematrix(:,1) == currentkiorder & formfunctionmonomialvariablematrix(:,2) == currentetaorder, 1))
                    continue;
                end
                
                if currentetaorder > 0
                    fprintf(fileid, ', ');
                end
                fprintf(fileid, '{');
                for currentphiorder = 0:maxphiorder
                    
                    if isempty(find(formfunctionmonomialvariablematrix(:,1) == currentkiorder & formfunctionmonomialvariablematrix(:,2) == currentetaorder & formfunctionmonomialvariablematrix(:,3) == currentphiorder, 1))
                        continue;
                    end
            
                    if currentphiorder > 0
                        fprintf(fileid, ', ');
                    end
                    
                    coefindex = find(formfunctionmonomialvariablematrix(:,1) == currentkiorder & formfunctionmonomialvariablematrix(:,2) == currentetaorder & formfunctionmonomialvariablematrix(:,3) == currentphiorder, 1);
                    currentcoefstring = char(vpa(  formfunctioncoefficientmatrix(i, coefindex)  , 20));
                    
                    fprintf(fileid, '%s', currentcoefstring);
                end
                fprintf(fileid, '}');
            end
            fprintf(fileid, '}');
        end
        fprintf(fileid, '});\n');    
    end
    fprintf(fileid, '               break;\n');
    
end


fprintf(fileid, '          default:\n');

fprintf(fileid, '               std::cout << "Error in ''%s'' namespace: Lagrange form functions of order %d and above not defined" << std::endl;\n', filename, order+1);
fprintf(fileid, '               abort();\n');
fprintf(fileid, '               break;\n');




fprintf(fileid, '     }\n');
fprintf(fileid, '\n');    
fprintf(fileid, '     return formfunctionpoly;\n');    
fprintf(fileid, '}');

fclose(fileid);
     
end
