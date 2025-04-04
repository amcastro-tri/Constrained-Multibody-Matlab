function N = kin_map(rpy)   
    % Diebel, Eq. 76
    theta = rpy(2); psi = rpy(3);
    sth = sin(theta); cth = cos(theta);
    sps = sin(psi); cps = cos(psi);
    N = [     cps,     sps, 0;
         -cth*sps, cth*cps, 0;
          cps*sth, sps*sth, cth]/cth;
end