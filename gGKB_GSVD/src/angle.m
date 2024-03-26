function [theta] = angle(v1,v2)
% Compute the sine of the angle of two vector,
% avoid computing sqrt(1-cos^2), since the loss of accuracy due 
% to root square when cos is approximate to 1.

[m1,n1] = size(v1);  [m2,n2] = size(v2);

if(n1 ~= 1 || n2 ~= 1 || m1 ~= m2)
    return;
end

v1 = v1/norm(v1);  v2 = v2/norm(v2);

e_proj = (eye(m1)-v1*v1')*v2;
theta = norm(e_proj);

end