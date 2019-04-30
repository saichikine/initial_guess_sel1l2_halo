function [zero_val,isterminal,direction] = poincare_event(t,X, plane, dist, L_point_num)

if strcmpi(plane,"xy")
    zero_val = X(3)-dist;
elseif strcmpi(plane,"yz")
    if L_point_num==2
        zero_val = X(1)-dist;
    elseif L_point_num ==1
        zero_val = dist-X(1);
    end
elseif strcmpi(plane,"xz")
    zero_val = X(2)-dist;
else
    error("this is bad")
end

isterminal = 1;  % Halt integration 
direction = -1;   % Decreasing (one sided poincare map)