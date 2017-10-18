function [ki, eta, phi] = getnodesinteriortotriangularface(nodecoord, order, takeall)

    % Default value:
    ki = [];
    eta = [];
    phi = [];
    
    % Number of nodes in each direction, including the exterior nodes:
    numnodes = (order+1);
    
    % In case there is no single node interior to the face:
    if numnodes <= 3 && takeall == 0
        return;
    end
    
    % 12 is the direction node1 -> node2, 23 the direction node2 -> node3, 31 the direction node3 to node 1:
    edge12length = nodecoord(:,2) - nodecoord(:,1);
    edge23length = nodecoord(:,3) - nodecoord(:,2);
    edge13length = nodecoord(:,3) - nodecoord(:,1);
    step12 = edge12length/(numnodes-1);
    step23 = edge23length/(numnodes-1);
    step13 = edge13length/(numnodes-1);
    
    % In case there is only a single central node interior to the face:
    if numnodes == 4 && takeall == 0
        centralnode = nodecoord(:,1) + step12 + step13;
        ki = centralnode(1);
        eta = centralnode(2);
        phi = centralnode(3);
        return;
    end

    % In case there is more than a single node interior to the face:
    % - redefine the 3 nodes as the ones one layer deeper in the face
    % - set order to order-3
    % - treat the corner and edge nodes
    % - do a recursive call
    
    if (takeall == 0)
        
        nodecoord(:,1) = nodecoord(:,1) + step12 + step13;
        nodecoord(:,2) = nodecoord(:,2) - step12 + step23;
        nodecoord(:,3) = nodecoord(:,3) - step13 - step23;

        order = order - 3;
        
    end
    
    % Take the corner nodes:
    ki = nodecoord(1,1:3);
    eta = nodecoord(2,1:3);
    phi = nodecoord(3,1:3);
    
    % Take the nodes interior to the edges:
    numedges = 3;
    edgedef = [1 2; 2 3; 3 1];
    for edge = 1:numedges
        [kiinneredge, etainneredge, phiinneredge] = getnodesinteriortoedge([ki(edgedef(edge,1)) ki(edgedef(edge,2)); eta(edgedef(edge,1)) eta(edgedef(edge,2)); phi(edgedef(edge,1)) phi(edgedef(edge,2))], order, 0);
        ki = [ki kiinneredge];
        eta = [eta etainneredge];
        phi = [phi phiinneredge];
    end
    
    % Perform the recursive call:
    [kirecursive, etarecursive, phirecursive] = getnodesinteriortotriangularface(nodecoord, order, 0);
    ki = [ki kirecursive];
    eta = [eta etarecursive];
    phi = [phi phirecursive];
    
end