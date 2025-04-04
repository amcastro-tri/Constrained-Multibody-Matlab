function S = skew(w)
    % Computes the skew matrix of a vector v, v× = skew(v), defined such
    % that skew(v)*a = v×a, the cross product with a.
    S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

