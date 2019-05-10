function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, y_k, proc_f, proc_Q, meas_h, meas_R)
    % PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
    %
    % Input:
    %   X_kmin1     [n x N] Particles for state x in time k-1
    %   W_kmin1     [1 x N] Weights for state x in time k-1
    %   y_k         [m x 1] Measurement vector for time k
    %   proc_f      Handle for process function f(x_k-1)
    %   proc_Q      [n x n] process noise covariance
    %   meas_h      Handle for measurement model function h(x_k)
    %   meas_R      [m x m] measurement noise covariance
    %
    % Output:
    %   X_k         [n x N] Particles for state x in time k
    %   W_k         [1 x N] Weights for state x in time k

    % calculate p(x(i)_k|x(i)_{k-1})
    X_k = mvnrnd( proc_f(X_kmin1)' ,proc_Q )';  
    
    % calculate p(y_k|x(i)_k) 
    Wy = mvnpdf(y_k',meas_h(X_k)',meas_R)';
    
    % compute weights
    W_k = W_kmin1 .* Wy;
    W_k = W_k / sum(W_k);
end