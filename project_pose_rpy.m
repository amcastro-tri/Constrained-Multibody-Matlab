function [x_proj, theta_proj] = project_pose_rpy(x, theta, params)
    max_iters = 10; tol = 1e-12;
    for i = 1:max_iters
        [Phi, J] = constraint_function(x, theta, params);
        if abs(Phi) < tol, break; end
        %J * J'
        delta = -J' * ((J * J') \ Phi);
        x = x + delta(1:3);
        theta = theta + delta(4:6);
    end
    x_proj = x;
    theta_proj = theta;
end
