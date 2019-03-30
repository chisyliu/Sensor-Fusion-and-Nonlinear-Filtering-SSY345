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
    % where E[x|y] is the mean of mix of multiple Gaussian densities
    %
    % In case of a mixture of one-dimensional distributions with weights wi, 
    % means mui and variances si2, the total mean and variance will be:
    % E[x|y] = sum_i { wi * μi } = w' * μi

    xHat = w'*mu;
end