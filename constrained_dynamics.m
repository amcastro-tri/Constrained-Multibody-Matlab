clear
close all

% Physical parameters
L = 0.1;           % Cube's size, meters.
radius = 0.1;      % Radius of the circle we constraint on, meters.
rho = 1000;        % Cube's density, kg/m³.
dissipation = 0;   % Translational dissipation, N⋅s/m.
g = 9.81;          % Acceleration of gravity, m/s².
sim_time = 5;             % Total simulated time.

% Time stepping parameters.
dt = 0.01;         % Time step, in seconds. 
the_v = 1.0;       % rpy method in translational velocity.
the_w = 1.0;       % rpy method in angular velocity.
c = 1e-10;         % Regularization.
num_steps = round(sim_time/dt);   % Number of time steps.

% Constraint parameters. 
params.p_BC = 0.5 * [-1; -1; 1] * L;  % Corner in body frame.
params.radius = radius;               % Circle radius.

% Cube
m = rho*L^3;       % Mass
Ixx = m/6*L^2;     % Inertia moments.
Iyy = Ixx; Izz = Ixx;

% Inertia tensor, about CoM Bcm, expressed in body frame B.
I_Bcm_B = diag([Ixx, Iyy, Izz]); 

% Initial configuration (must satisfies constraint)
x = [1.5; 0.5; -0.5] * L;
rpy = [0; 0; 0];  % RPY
[Phi, J] = constraint_function(x, rpy, params)

% Project onto constraint manifold.
[x, rpy] = project_pose_rpy(x, rpy, params);
[Phi, J] = constraint_function(x, rpy, params)
nk = length(Phi);

% Initial velocities.
v = [0; 0.1; 0];
omega = [0; 0; 0];

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

vid = VideoWriter('cube_animation.avi', 'Motion JPEG AVI');  % Create video object
vid.FrameRate = 30;  % Adjust as needed
open(vid);

% Simulation loop
for k = 1:num_steps
    t = (k-1) * dt;
    R = rpy2rotm(rpy);
    I_world = R * I_Bcm_B * R';  % world-frame inertia

    % External forces and torques
    f_ext = [0; 0; -m*g];
    tau_ext = [0; 0; 0];

    % Compute free motion velocities (no constraints).
    v_star =  (m * v + dt * f_ext) / (m + dt * dissipation);
    L_Bcm_W = I_world * omega;
    C_W = cross(omega, L_Bcm_W);
    omega_star = omega + dt * (I_world \ (tau_ext - C_W));

    % Solve for Lagrange multipliers.
    [Phi, J] = constraint_function(x, rpy, params);        
    M = blkdiag(m * eye(3), I_world); % Mass matrix
    Phi_star = Phi + dt * J * [v_star; omega_star];
    MiJT = M \ J';
    W = J * MiJT;
    Wreg = W + c/dt * eye(nk);
    lambda = -Wreg \ Phi_star / dt;

    % Forces and moments on the body.
    forces = J' * lambda;

    % Incorporate Lagrange multipliers in the full velocity update.
    delta_v = MiJT * lambda;
    v_next = v_star + delta_v(1:3);
    omega_next = omega_star + delta_v(4:6);

    % Integrate positions. We use the rpy method.
    x_next = x + dt * (the_v * v_next + (1-the_v) * v);
    rpy_dot = kin_map(rpy) * (the_w * omega_next + (1-the_w) * omega);
    rpy_next = rpy + dt * rpy_dot;

    % Project to constraint manifold
    [x_next, rpy_next] = project_pose_rpy(x_next, rpy_next, params);

    % Update graphics
    R_next = rpy2rotm(rpy_next);
    C_w = R_next * C_b + x_next;
    set(cube_patch, 'Vertices', C_w');
    drawnow;
    axis([-2 3 -2 2 -2 2]*radius)
    
    frame = getframe(gcf);  % Capture current figure
    writeVideo(vid, frame);   % Write to file

    % Update state
    x = x_next;
    rpy = rpy_next;
    v = v_next;
    omega = omega_next;

    T(k) = k *dt;
    K(k) = 0.5 * omega' * I_world * omega + 0.5 * m * v'*v;
    V(k) = m * g * x(3);
    E(k) = K(k) + V(k);
    U(:,k) = rpy;
    F(:,k) = forces;

    pause(dt/10);
end

close(vid);  % Finalize the MP4 file


