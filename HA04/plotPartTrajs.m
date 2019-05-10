
function plotPartTrajs(k, Xk, Xkmin1, ~, j)
%PLOTPARTTRAJS Summary of this function goes here
%   Plots lines between ith sample of Xk and j(i)th sample of Xk-1. When 
%   repeated during a particle filter execution, this will produce particle 
%   trajectories illustration over time.
%
%   This function is intended to be passed as a function handle into your
%   particle filter function.
%
% Inputs:
%   k           time instance index
%   Xk          [n x N] N particles of dimension n to approximate p(x_k).
%   Xkmin1      [n x N] N particles of dimension n to approximate p(x_k-1).
%   Wk          [1 x N] Corresponding weights.
%   j           Index vector such that Xk(:,i) = Xkmin1(:,j(i))

    if (size(Xk,2) <= 50) % At most 50 particles may be plotted
        for i = 1:size(Xk,2) % loop through all particles
            plot([k-1 k], [Xkmin1(1,j(i)) Xk(1,i)]);
            hold on 
        end
        title(['Particle trajectories up to time k=', num2str(k)]);
        pause(0.05);
    else
        disp('Too many particles to plot!'); 
    end
end


