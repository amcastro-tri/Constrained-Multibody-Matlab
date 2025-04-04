function [x_proj, rpy_proj] = project_pose_rpy(x, rpy, params)
    max_iters = 10; tol = 1e-12;
    for i = 1:max_iters
        [Phi, J] = constraint_function(x, rpy, params);
        if abs(Phi) < tol, break; end
        %J * J'
        delta = -J' * ((J * J') \ Phi);
        x = x + delta(1:3);
        rpy = rpy + delta(4:6);
    end
    x_proj = x;
    rpy_proj = rpy;
end
