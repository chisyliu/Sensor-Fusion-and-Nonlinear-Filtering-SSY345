function [xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, Q, h, R, sigmaPoints, type)
    % NONLINRTSSMOOTHER Filters measurement sequence Y using a 
    % non-linear Kalman filter. 
    %
    % Input:
    %   Y           [m x N] Measurement sequence for times 1,...,N
    %   x_0         [n x 1] Prior mean for time 0
    %   P_0         [n x n] Prior covariance
    %   f                   Motion model function handle
    %   T                   Sampling time
    %   Q           [n x n] Process noise covariance
    %   S           [n x N] Sensor position vector sequence
    %   h                   Measurement model function handle
    %   R           [n x n] Measurement noise covariance
    %   sigmaPoints Handle to function that generates sigma points.
    %   type        String that specifies type of non-linear filter/smoother
    %
    % Output:
    %   xf          [n x N]     Filtered estimates for times 1,...,N
    %   Pf          [n x n x N] Filter error convariance
    %   xp          [n x N]     Predicted estimates for times 1,...,N
    %   Pp          [n x n x N] Filter error convariance
    %   xs          [n x N]     Smoothed estimates for times 1,...,N
    %   Ps          [n x n x N] Smoothing error convariance

    
    %% Parameters
    N = size(Y,2);
    n = length(x_0);
    m = size(Y,1);
    
    
    %% Forward filtering
    
%     f2 = @(x) f(x,T);
%     h2 = @(x) h(x,S(:,1));
    [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type);
    
    %% Backward filtering

    % initialize outputs
    xs(:,N) = xf(:,N);
    Ps(:,:,N) = Pf(:,:,N);

    % backward recursion - last filter output is already smoothed
    for k=N-1:-1:1
        [xs(:,k), Ps(:,:,k)] = nonLinRTSSupdate(xs(:,k+1), Ps(:,:,k+1), xf(:,k), Pf(:,:,k), xp(:,k+1),  Pp(:,:,k+1), f, sigmaPoints, type);
    end
    
end
