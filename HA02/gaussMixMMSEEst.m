function [ xHat ] = gaussMixMMSEEst( w, mu, sigma2 )
    %GAUSSMIXMMSEEST calculates the MMSE estimate from a Gaussian mixture
    %density with multiple components.
    %
    %Input
    %   W           Vector of all the weights
    %   MU          Vector containing the means of all components
    %   SIGMA2      Vector containing the variances of all components
    %
    %Output
    %   xHat        MMSE estimate

    % The MMSE estimator of x, given some observation y, is equal to E[x|y]
    % where E[x|y] is the mean of mix of multiple Gaussian densities p(x|y):
    % 
    % p(x|y) = w1 * N1(x;mu,var) +...+ wn * Nn(x;mu,var),   |w|_1=1, 0<wi<1
    %
    % In case of a mixture of one-dimensional Gaussian distributions
    % weighted by wi, with means mui and variances si2, the MMSE will be:
    %
    % x_{MMSE} = E[x|y] = sum_i { wi * μi } = w' * μi

    xHat = w(:)'*mu(:);
end