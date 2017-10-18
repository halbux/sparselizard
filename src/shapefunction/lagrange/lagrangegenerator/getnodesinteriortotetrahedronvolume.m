function [ki, eta, phi] = getnodesinteriortotetrahedronvolume(nodecoord, order, takeall)

    % Default value:
    ki = [];
    eta = [];
    phi = [];
    
    % Number of nodes in each direction, including the exterior nodes:
    numnodes = (order+1);
    
    % In case there is no single node interior to the face:
    if numnodes <= 4 && takeall == 0
        return;
    end
    
    % 12 is the direction node1 -> node2, ...:
    edge12length = nodecoord(:,2) - nodecoord(:,1);
    edge23length = nodecoord(:,3) - nodecoord(:,2);
    edge13length = nodecoord(:,3) - nodecoord(:,1);
    edge14length = nodecoord(:,4) - nodecoord(:,1);
    edge24length = nodecoord(:,4) - nodecoord(:,2);
    edge34length = nodecoord(:,4) - nodecoord(:,3);
    step12 = edge12length/(numnodes-1);
    step23 = edge23length/(numnodes-1);
    step13 = edge13length/(numnodes-1);
    step14 = edge14length/(numnodes-1);
    step24 = edge24length/(numnodes-1);
    step34 = edge34length/(numnodes-1);
    
    % In case there is only a single central node interior to the face:
    if numnodes == 5 && takeall == 0
        centralnode = nodecoord(:,1) + step12 + step13 + step14;
        ki = centralnode(1);
        eta = centralnode(2);
        phi = centralnode(3);
        return;
    end

    % In case there is more than a single node interior to the face:
    % - redefine the 4 nodes as the ones one layer deeper in the face
    % - set order to order-3
    % - treat the corner and edge nodes
    % - do a recursive call
    
    if (takeall == 0)
    
        nodecoord(:,1) = nodecoord(:,1) + step12 + step13 + step14;
        nodecoord(:,2) = nodecoord(:,2) - step12 + step23 + step24;
        nodecoord(:,3) = nodecoord(:,3) - step13 - step23 + step34;
        nodecoord(:,4) = nodecoord(:,4) - step14 - step24 - step34;

        order = order - 4;
        
    end
    
    % Take the corner nodes:
    ki = nodecoord(1,1:4);
    eta = nodecoord(2,1:4);
    phi = nodecoord(3,1:4);
    
    % Treat the nodes interior to the edges:
    numedges = 6;
    edgedef = [1 2; 2 3; 3 1; 4 1; 4 3; 4 2];
    for edge = 1:numedges
        [kiinneredge, etainneredge, phiinneredge] = getnodesinteriortoedge([ki(edgedef(edge,1)) ki(edgedef(edge,2)); eta(edgedef(edge,1)) eta(edgedef(edge,2)); phi(edgedef(edge,1)) phi(edgedef(edge,2))], order, 0);
        ki = [ki kiinneredge];
        eta = [eta etainneredge];
        phi = [phi phiinneredge];
    end
    
    % Take the nodes interior to the faces:
    numfaces = 4;
    facedef = [1 3 2; 1 2 4; 1 4 3; 4 2 3];
    for face = 1:numfaces
        [kiinnerface, etainnerface, phiinnerface] = getnodesinteriortotriangularface([ki(facedef(face,1)) ki(facedef(face,2)) ki(facedef(face,3)); eta(facedef(face,1)) eta(facedef(face,2)) eta(facedef(face,3)); phi(facedef(face,1)) phi(facedef(face,2)) phi(facedef(face,3))], order, 0);
        ki = [ki kiinnerface];
        eta = [eta etainnerface];
        phi = [phi phiinnerface];
    end
    
    % Perform the recursive call:
    [kirecursive, etarecursive, phirecursive] = getnodesinteriortotetrahedronvolume(nodecoord, order, 0);
    ki = [ki kirecursive];
    eta = [eta etarecursive];
    phi = [phi phirecursive];
    
end