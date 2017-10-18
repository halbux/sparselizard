function [ki, eta, phi] = getnodesinteriortoprismvolume(nodecoord, order, takeall)

    % Default value:
    ki = [];
    eta = [];
    phi = [];
    
    % Take the corner nodes:
    ki = nodecoord(1,1:6);
    eta = nodecoord(2,1:6);
    phi = nodecoord(3,1:6);
    
    % Treat the nodes interior to the edges:
    numedges = 9;
    edgedef = [1 2; 1 3; 1 4; 2 3; 2 5; 3 6; 4 5; 4 6; 5 6];
    for edge = 1:numedges
        [kiinneredge, etainneredge, phiinneredge] = getnodesinteriortoedge([ki(edgedef(edge,1)) ki(edgedef(edge,2)); eta(edgedef(edge,1)) eta(edgedef(edge,2)); phi(edgedef(edge,1)) phi(edgedef(edge,2))], order, 0);
        ki = [ki kiinneredge];
        eta = [eta etainneredge];
        phi = [phi phiinneredge];
    end
    
    % Take the nodes interior to the faces:
    numtriangularfaces = 2;
    facedef = [1 3 2; 4 5 6];
    for face = 1:numtriangularfaces
        [kiinnerface, etainnerface, phiinnerface] = getnodesinteriortotriangularface([ki(facedef(face,1)) ki(facedef(face,2)) ki(facedef(face,3)); eta(facedef(face,1)) eta(facedef(face,2)) eta(facedef(face,3)); phi(facedef(face,1)) phi(facedef(face,2)) phi(facedef(face,3))], order, 0);
        ki = [ki kiinnerface];
        eta = [eta etainnerface];
        phi = [phi phiinnerface];
    end
    numtquadrangularfaces = 3;
    facedef = [1 2 5 4; 1 4 6 3; 2 3 6 5];
    for face = 1:numtquadrangularfaces
        [kiinnerface, etainnerface, phiinnerface] = getnodesinteriortoquadrangularface([ki(facedef(face,1)) ki(facedef(face,2)) ki(facedef(face,3)) ki(facedef(face,4)); eta(facedef(face,1)) eta(facedef(face,2)) eta(facedef(face,3)) eta(facedef(face,4)); phi(facedef(face,1)) phi(facedef(face,2)) phi(facedef(face,3)) phi(facedef(face,4))], order, 0);
        ki = [ki kiinnerface];
        eta = [eta etainnerface];
        phi = [phi phiinnerface];
    end
    
    % Treat the nodes interior to the volume:
    
    % Number of nodes in each direction, including the exterior nodes:
    numnodes = (order+1);
    
    % In case there is no single node interior to the volume:
    if numnodes <= 3 && takeall == 0
        return;
    end
    
    % Get all nodes interior to face [1 2 3] (don't known why not 132 but that's the way it is in gmsh...)
    [kiinnerface123, etainnerface123, phiinnerface123] = getnodesinteriortotriangularface([ki(1) ki(2) ki(3); eta(1) eta(2) eta(3); phi(1) phi(2) phi(3)], order, 0);

    edge14length = nodecoord(:,4) - nodecoord(:,1);
    step14 = edge14length/(numnodes-1);
    
    for i = 1:length(kiinnerface123)
        
        % The face right above the 123 face:
        [kiinnervolume, etainnervolume, phiinnervolume] = getnodesinteriortoedge([kiinnerface123(i) kiinnerface123(i); etainnerface123(i) etainnerface123(i); phiinnerface123(i)+step14(3) phiinnerface123(i)+edge14length(3)-step14(3)], order-2, 1);

        ki = [ki kiinnervolume];
        eta = [eta etainnervolume];
        phi = [phi phiinnervolume];
        
    end
end