function [mu, sigma2] = posteriorGaussian(mu_x, sigma2_x, y, sigma2_r)
    % Calculates the posterior p(x|y) which is proportional to p(x,y)
    %
    % posteriorGaussian performs a single scalar measurement update with a
    % measurement model which is simply "y = x + noise".
    %
    % Input
    %   MU_X            The mean of the (Gaussian) prior density.
    %   SIGMA2_X        The variance of the (Gaussian) prior density.
    %   SIGMA2_R        The variance of the measurement noise.
    %   Y               The given measurement.
    %
    % Output
    %   MU              The mean of the (Gaussian) posterior distribution
    %   SIGMA2          The variance of the (Gaussian) posterior distribution

    % % Calculate Joint distribution p(x,y) which is proportional to p(x|y)
    % syms mux sx sr x y
    % [mu_xy, Q_xy] = jointGaussian(mux, sx, sr);

    % % Express Gaussian distribution p(x,y)= \propto p(x|y) in terms of x (given y)
    % syms munew snew c
    % eqq1 = ([x;y] - mu_xy ).' * Q_xy^-1  * ([x;y] - mu_xy )     ;
    % eqq2 = (  x   - munew ).  * snew^-1  * (  x   - munew )  + c;
    % sol = solve( coeffs( eqq1 - eqq2, x) , [munew snew c]);

    % % Show results
    % simplify(sol.munew);  % -> (mux*sr + sx*y)/(sr + sx)
    % simplify(sol.snew);   % -> (sr*sx)/(sr + sx)
    % simplify(sol.c);      % -> (mux - y)^2/(sr + sx)

    % apply results from symbolic calculations
    mu  = (mu_x*sigma2_r + sigma2_x*y)/(sigma2_r + sigma2_x);
    sigma2 = (1/sigma2_r + 1/sigma2_x)^-1;
end