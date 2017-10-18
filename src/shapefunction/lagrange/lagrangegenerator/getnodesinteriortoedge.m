function [ki, eta, phi] = getnodesinteriortoedge(nodecoord, order, takeall)

    % Default value:
    ki = [];
    eta = [];
    phi = [];
    
    if (takeall == 1)
        % Take the corner nodes:
        ki = nodecoord(1,1:2);
        eta = nodecoord(2,1:2);
        phi = nodecoord(3,1:2);
    end

    % Number of nodes on the edge, including the corner nodes:
    numnodes = (order+1);
    
    % In case there is no single node interior to the edge:
    if numnodes == 2
        return;
    end
    
    edgelength = nodecoord(:,2) - nodecoord(:,1);
    step = edgelength/(numnodes-1);

    for i = 1:numnodes - 2
        ki = [ki nodecoord(1,1) + i*step(1)];
        eta = [eta nodecoord(2,1) + i*step(2)];
        phi = [phi nodecoord(3,1) + i*step(3)];
    end
    
end




