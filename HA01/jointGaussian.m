function [mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)
    % JointGaussian calculates the joint Gaussian density p([y;x])
    %
    % y = x + r
    %               x ~ N(mu_x,sigma2_x)
    %               r ~ N(0   ,sigma2_r)
    %
    %  =>  p([y;x]) = p(A*[x;r] + b) = p( [1 0; 1 1]*[x;r] + [0;0] )
    %
    % Input
    %   MU_X        Expected value of x
    %   SIGMA2_X    Covariance of x
    %   SIGMA2_R    Covariance of the noise r
    %
    % Output
    %   MU          Mean of joint density 
    %   SIGMA       Covariance of joint density

    % define linear transformations matrices [x;y] = A*[x;r] + b according to
    % problem 1.3a
    A_xr2xy  = [1 0; 1 1];
    b_xr2xy  = [0;0];

    % define mean of [x;r] = [E[x] ; E[r]] = [mu_x, mu_r=0]
    mu_xr    = [mu_x; 0];

    % define covariance of [x;r]. Since they are independent, 
    % then cov[x;r] = [ cov[x] 0; 0; cov[r] ] = diag(cov[x], cov[r])
    Sigma_xr = blkdiag(sigma2_x,sigma2_r);

    % calculate mean and cov of the new vector [x;y], which is obtained from a
    % linear transformation [x;y] = A*[x;r] + b
    [mu, Sigma] = affineGaussianTransform(mu_xr, Sigma_xr, A_xr2xy, b_xr2xy);
end