clear
close all

%function constrained_dynamics()
    % Parameters
    L = 0.1;
    radius = 0.1;
    rho = 1000;    
    dissipation = 0;
    g = 9.81;
%    I_body = [0.2 0.05 0.0;
%              0.05 0.25 0.01;
%              0.0 0.01 0.3];  % arbitrary symmetric inertia tensor

    % Cube
    m = rho*L^3;
    Ixx = m/6*L^2;
    Iyy = Ixx; Izz = Ixx;
    I_body = diag([Ixx, Iyy, Izz]);

    dt = 0.01;
    the_v = 1.0; the_w = 1.0;
    T = 5;
    c = 1e-10;  % Regularization.
    N = round(T/dt);

    % Initial state (satisfies constraint)
    x = [1.5; 0.5; -0.5] * L;
    theta = [0; 0; 0];  % RPY
    v = [0; 0.1; 0];
    %omega = 2 * pi * [0; 1; 1]/sqrt(2);
    omega = [0; 0; 0];

    % Constrained corner in body frame
    corner_b = 0.5 * [-1; -1; 1] * L;

    params.p_BC = corner_b;
    params.radius = radius;

    [Phi, J_full] = constraint_function(x, theta, params)
    nk = length(Phi);    

    % Cube definition
    C_b = L * 0.5 * [
        -1 -1 -1;
         1 -1 -1;
         1  1 -1;
        -1  1 -1;
        -1 -1  1;
         1 -1  1;
         1  1  1;
        -1  1  1
    ]';

    % Graphics
    figure; axis equal; grid on; hold on;
    xlabel('X'); ylabel('Y'); zlabel('Z'); view(3);
    th = linspace(0, 2*pi, 100);
    plot3(radius*cos(th), radius*sin(th), zeros(size(th)), 'k--', 'LineWidth', 1.2);
    cube_patch = patch('Faces', cube_faces(), 'Vertices', zeros(8,3), ...
        'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5);

    % Simulation loop
    for k = 1:N
        t = (k-1) * dt;
        R = rpy2rotm(theta);
        I_world = R * I_body * R';  % world-frame inertia

        % External forces and torques
        f_ext = [0; 0; -m*g];
        tau_ext = [0; 0; 0];

        % Tentative updates
        v_star =  (m * v + dt * f_ext) / (m + dt * dissipation);
        L_Bcm_W = I_world * omega;
        C_W = cross(omega, L_Bcm_W);
        omega_star = omega + dt * (I_world \ (tau_ext - C_W));

        % Constraint function and Jacobian
        [Phi, J_full] = constraint_function(x, theta, params);
        %fprintf("t: %g, phi: %g\n", t, norm(Phi));
        %gamma = -Phi / dt;

        % Mass matrix
        M = blkdiag(m * eye(3), I_world);

        % Constraint impulse
        Phi_star = Phi + dt * J_full * [v_star; omega_star];
        MiJT = M \ J_full';
        W = J_full * MiJT;
        Wreg = W + c/dt * eye(nk);
        lambda = -Wreg \ Phi_star / dt;
        delta_v = MiJT * lambda;
        %delta_v = zeros(6, 1);

        % Corrected velocities
        v_next = v_star + delta_v(1:3);
        omega_next = omega_star + delta_v(4:6);

        % Integrate
        x_next = x + dt * (the_v * v_next + (1-the_v) * v);
        theta_dot = T_rpy(theta) * (the_w * omega_next + (1-the_w) * omega);
        theta_next = theta + dt * theta_dot;

        % Project to constraint manifold
        [x_next, theta_next] = project_pose_rpy(x_next, theta_next, params);        

        % Update graphics
        R_next = rpy2rotm(theta_next);
        C_w = R_next * C_b + x_next;
        %C_w = C_b + x;
        set(cube_patch, 'Vertices', C_w');
        drawnow;
        axis([-2 3 -2 2 -2 2]*radius)
        pause(dt/10);

        % Update state
        x = x_next;
        theta = theta_next;
        v = v_next;
        omega = omega_next;

        T(k) = k *dt;    
        K(k) = 0.5 * omega' * I_world * omega + 0.5 * m * v'*v;
        V(k) = m * g * x(3);
        E(k) = K(k) + V(k);
        U(:,k) = theta;

    end


