function R = rpy2rotm(rpy)
    % Computes the rotation matrix R_WB for the pose of body B in the world
    % frame W, given rpy angles.
    phi = rpy(1); the = rpy(2); psi = rpy(3);
    Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
    Ry = [cos(the), 0, sin(the); 0, 1, 0; -sin(the), 0, cos(the)];
    Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
    R = Rz * Ry * Rx;
end
