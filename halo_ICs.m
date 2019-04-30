function ICs = halo_computeplot(Az, mu, L, L_point_string, halo_pole)
    % Returns initial conditions for Halo orbit for a mu value and fixed Az
    
    %% Input checking
    
    if ~isstring(halo_pole)
        error('Fifth argument must be a string.')
    end
    
    if (~strcmp(lower(halo_pole),'south') && ~strcmp(lower(halo_pole),'north'))
        error('Fifth argument must be either "south" or "north"')
    end
    
    if ~isstring(L_point_string)
        error('Fourth argument must be a string.')
    end
    
    if (~strcmp(lower(L_point_string),'l1') && ~strcmp(lower(L_point_string),'l2'))
        error('Fourth argument must be either "L1" or "L2"')
    end
    
    %% L1 or L2?
    
    L_points = lagrangePoints(mu)
    
    if strcmp(lower(L_point_string),'l1')
        L_point = L_points(:,1);
    elseif strcmp(lower(L_point_string),'l2')
        L_point = L_points(:,2);
    else
        error('You should not be here.')
    end
    
    
end