function [ki, eta, phi] = getnodesinteriortoquadrangularface(nodecoord, order, takeall)

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
    
    % 12 is the direction node1 -> node2, 23 the direction node2 -> node3:
    edge12length = nodecoord(:,2) - nodecoord(:,1);
    edge23length = nodecoord(:,3) - nodecoord(:,2);
    step12 = edge12length/(numnodes-1);
    step23 = edge23length/(numnodes-1);
    
    % In case there is only a single central node interior to the face:
    if numnodes == 3 && takeall == 0
        centralnode = nodecoord(:,1) + step12 + step23;
        ki = centralnode(1);
        eta = centralnode(2);
        phi = centralnode(3);
        return;
    end

    % In case there is more than a single node interior to the face:
    % - redefine the 4 nodes as the ones one layer deeper in the face
    % - set order to order-2
    % - treat the corner and edge nodes
    % - do a recursive call
    
    if (takeall == 0)
        
        nodecoord(:,1) = nodecoord(:,1) + step12 + step23;
        nodecoord(:,2) = nodecoord(:,2) - step12 + step23;
        nodecoord(:,3) = nodecoord(:,3) - step12 - step23;
        nodecoord(:,4) = nodecoord(:,4) + step12 - step23;

        order = order - 2;
        
    end
    
    % Take the corner nodes:
    ki = nodecoord(1,1:4);
    eta = nodecoord(2,1:4);
    phi = nodecoord(3,1:4);
    
    % Take the nodes interior to the edges:
    numedges = 4;
    edgedef = [1 2; 2 3; 3 4; 4 1];
    for edge = 1:numedges
        [kiinneredge, etainneredge, phiinneredge] = getnodesinteriortoedge([ki(edgedef(edge,1)) ki(edgedef(edge,2)); eta(edgedef(edge,1)) eta(edgedef(edge,2)); phi(edgedef(edge,1)) phi(edgedef(edge,2))], order, 0);
        ki = [ki kiinneredge];
        eta = [eta etainneredge];
        phi = [phi phiinneredge];
    end
    
    % Perform the recursive call:
    [kirecursive, etarecursive, phirecursive] = getnodesinteriortoquadrangularface(nodecoord, order, 0);
    ki = [ki kirecursive];
    eta = [eta etarecursive];
    phi = [phi phirecursive];
    
end