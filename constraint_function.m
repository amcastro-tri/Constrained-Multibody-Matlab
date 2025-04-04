function [Phi, J] = constraint_function(x, theta, params)
    % Body-fixed constrained point
    %c_b = 0.5 * [-1; -1; 1];
    c_b = params.p_BC;
    r = params.radius;
    r2 = r * r;
    
    R = rpy2rotm(theta);
    c_w = x + R * c_b;

    % Constraint: point lies on unit circle in XY
    Phi = zeros(2, 1);
    Phi(1) = c_w(1)^2 + c_w(2)^2 - r2;
    Phi(2) = c_w(3);

    % Derivative of constraint wrt c_w
    dPhi_dc = [2*c_w(1), 2*c_w(2), 0;
                      0,        0, 1];

    % Derivatives wrt generalized velocities
    Jx = eye(3);
    Jr = -skew(R * c_b);
    J = dPhi_dc * [Jx, Jr];  % 2x6 Jacobian
end
