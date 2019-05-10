function plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, ax)
%PLOTPOSTPDF Plots blurred pdf for a PF posterior, and plots a Kalman
% posterior to compare with.
%
%   This function is intended to be used as a function handle for a
%   compatible particle filter function. It is meant to be called each time
%   the particle filter has updated the particles (but before any
%   resampling has been carried out.)
%
%   To use it in your filter you should first compute xf, Pf, and set
%   bResample, sigma and ax.
%   Then define a function handle
%       plotFunc_handle  = @(k, Xk, Xkmin1, Wk, j) ...
%                          (plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, ax))
%   Then call your PF-function with plotFunc_handle as plotFunc argument.
%
% Inputs:
%   k           time instance index
%   Xk          [n x N] N particles of dimension n to approximate p(x_k).
%   Wk          [1 x N] Corresponding weights.
%   xf          [n x K] Filter posteriors for some filter to compare with
%   Pf          [n x n x K] Filter posterior covariances for ^
%   bResample   Flag for resampling. False: do not resample, True: resample
%   sigma       Controls the kernel width for blurring particles.
%   ax          [xmin xmax ymin ymax] Used for setting the x-axis limits to
%               a value that doesn't change through iterations of the PF
%               filter.

    N = size(Xk,2);

    % Let us first determine the x-interval of interest:
    xmin =    min(Xk(1,:)); %ax(1);
    xmax =    max(Xk(1,:)); %ax(2); 
    X    =    linspace(xmin-(xmax-xmin)/3, xmax+(xmax-xmin)/3, 800);

    % We can now construct a continuous approximation to the posterior
    % density by placing a Gaussian kernel around each particle
    pApprox = zeros(size(X));   % A vector that will contain the pdf values

    if bResample
        sigma=(xmax-xmin)/sqrt(N);
    end
    
    for i = 1 : N
        pApprox = pApprox + Wk(1,i)*normpdf(Xk(1,i), X, sigma);
    end

    % We are now ready to plot the densities
    
    % figure;
    set(gcf, 'Name', ['p_',num2str(k), '_', 'SIR']);
    % clf
    
    plot(X, pApprox, 'LineWidth', 2)   % This is the PF approximation
    hold on
    plot(X, normpdf(xf(1,k), X, sqrt(Pf(1,1,k))), 'r-.', 'LineWidth', 2) % KF posterior density
    legend('Particle filter approximation', 'Kalman filter', 'Location', 'southwest')
    title(['p(x_k |  y_{1:k}), k=', num2str(k)])
    hold off;
    pause()
end


