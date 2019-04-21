function [mu_y, Sigma_y] = affineGaussianTransform(mu_x, Sigma_x, A, b)
    %affineTransformGauss calculates the mean and covariance of y, the 
    %transformed variable, exactly when the function, f, is defined as 
    %y = f(x) = Ax + b, where A is a matrix, b is a vector of the same 
    %dimensions as y, and x is a Gaussian random variable.
    %
    %Input
    %   MU_X        [n x 1] Expected value of x.
    %   SIGMA_X     [n x n] Covariance of x.
    %   A           [m x n] Linear transform matrix.
    %   B           [m x 1] Constant part of the affine transformation.
    %
    %Output
    %   MU_Y        [m x 1] Expected value of y.
    %   SIGMA_Y     [m x m] Covariance of y.


    % E[A*x + b] = A*E[x] + b = A*mu_x + b
    mu_y = A * mu_x + b;
    % Cov[A*x + b] = A*Cov[x]*A^T
    Sigma_y = A * Sigma_x * A';
end