function plotPartTrajs(k, Xk, Xkmin1, Wk, j, bAlpha)
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
    
    hold on;
    N = size(Xk,2);
    if ( N <= 50) % At most 50 particles may be plotted
        for i = 1:N % loop through all particles
            if bAlpha
                alph = min(1,max(0, Wk(i)*N ));
            else
                alph = 0.2;
            end
            plot([k-1 k], [Xkmin1(1,j(i)) Xk(1,i)], 'LineWidth',1,'Color',[0 0 0 alph]);
        end
        title(['Particle trajectories up to time k=', num2str(k)]);
        drawnow();
    else
        disp('Too many particles to plot!'); 
    end
end


