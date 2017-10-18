function [ki, eta, phi] = getnodesinteriortohexahedronvolume(nodecoord, order, takeall)

    % Default value:
    ki = [];
    eta = [];
    phi = [];
    
    % Number of nodes in each direction, including the exterior nodes:
    numnodes = (order+1);
    
    % In case there is no single node interior to the face:
    if numnodes == 2 && takeall == 0
        return;
    end
    
    % 12 is the direction node1 -> node2, 23 the direction node2 -> node3, 15 the direction node1 -> node5:
    edge12length = nodecoord(:,2) - nodecoord(:,1);
    edge23length = nodecoord(:,3) - nodecoord(:,2);
    edge15length = nodecoord(:,5) - nodecoord(:,1);
    step12 = edge12length/(numnodes-1);
    step23 = edge23length/(numnodes-1);
    step15 = edge15length/(numnodes-1);
    
    % In case there is only a single central node interior to the face:
    if numnodes == 3 && takeall == 0
        
        centralnode = nodecoord(:,1) + step12 + step23 + step15;
        ki = centralnode(1);
        eta = centralnode(2);
        phi = centralnode(3);
        return;
    end

    % In case there is more than a single node interior to the face:
    % - redefine the 8 nodes as the ones one layer deeper in the face
    % - set order to order-2
    % - treat the corner and edge nodes
    % - do a recursive call
    
    if (takeall == 0)
  
        nodecoord(:,1) = nodecoord(:,1) + step12 + step23 + step15;
        nodecoord(:,2) = nodecoord(:,2) - step12 + step23 + step15;
        nodecoord(:,3) = nodecoord(:,3) - step12 - step23 + step15;
        nodecoord(:,4) = nodecoord(:,4) + step12 - step23 + step15;
        nodecoord(:,5) = nodecoord(:,5) + step12 + step23 - step15;
        nodecoord(:,6) = nodecoord(:,6) - step12 + step23 - step15;
        nodecoord(:,7) = nodecoord(:,7) - step12 - step23 - step15;
        nodecoord(:,8) = nodecoord(:,8) + step12 - step23 - step15;

        order = order - 2;
        
    end
    
    % Take the corner nodes:
    ki = nodecoord(1,1:8);
    eta = nodecoord(2,1:8);
    phi = nodecoord(3,1:8);
    
    % Treat the nodes interior to the edges:
    numedges = 12;
    
    edgedef = [1 2; 1 4; 1 5; 2 3; 2 6; 3 4; 3 7; 4 8; 5 6; 5 8; 6 7; 7 8];
    for edge = 1:numedges
        [kiinneredge, etainneredge, phiinneredge] = getnodesinteriortoedge([ki(edgedef(edge,1)) ki(edgedef(edge,2)); eta(edgedef(edge,1)) eta(edgedef(edge,2)); phi(edgedef(edge,1)) phi(edgedef(edge,2))], order, 0);
        ki = [ki kiinneredge];
        eta = [eta etainneredge];
        phi = [phi phiinneredge];
    end
    
    % Take the nodes interior to the faces:
    numfaces = 6;

    facedef = [1 4 3 2; 1 2 6 5; 1 5 8 4; 2 3 7 6; 3 4 8 7; 5 6 7 8];
    for face = 1:numfaces
        [kiinnerface, etainnerface, phiinnerface] = getnodesinteriortoquadrangularface([ki(facedef(face,1)) ki(facedef(face,2)) ki(facedef(face,3)) ki(facedef(face,4)); eta(facedef(face,1)) eta(facedef(face,2)) eta(facedef(face,3)) eta(facedef(face,4)); phi(facedef(face,1)) phi(facedef(face,2)) phi(facedef(face,3)) phi(facedef(face,4))], order, 0);
        ki = [ki kiinnerface];
        eta = [eta etainnerface];
        phi = [phi phiinnerface];
    end
    
    % Perform the recursive call:
    [kirecursive, etarecursive, phirecursive] = getnodesinteriortohexahedronvolume(nodecoord, order, 0);
    ki = [ki kirecursive];
    eta = [eta etarecursive];
    phi = [phi phirecursive];
    
end