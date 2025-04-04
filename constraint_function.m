function [Phi, J] = constraint_function(p_WB, rpy, params)    
    p_BC = params.p_BC;
    r = params.radius;
    r2 = r * r;
    
    R_WB = rpy2rotm(rpy);
    p_BC_W = R_WB * p_BC;
    p_WC = p_WB + p_BC_W;

    % Constraint: point lies on unit circle in XY, at Z = 0.
    % This is two constraint equations.
    Phi = zeros(2, 1);
    Phi(1) = p_WC(1)^2 + p_WC(2)^2 - r2;
    Phi(2) = p_WC(3);

    % Derivative of constraint wrt p_WC
    dPhi_dp = [2*p_WC(1), 2*p_WC(2), 0;
                      0,        0, 1];
    
    % Corner velocity Jacobian, i.e. the corner velocity in the world is:
    %   v_WC = J_WC ⋅ v.
    Jx = eye(3);
    Jw = -skew(p_BC_W);
    J_WC = [Jx, Jw];

    % Geometric Jacobian.
    % i.e dPhi/dt = J⋅v + ∂Φ/∂t, with v generalized velocities.
    J = dPhi_dp * J_WC;  % 2x6 Jacobian
end
