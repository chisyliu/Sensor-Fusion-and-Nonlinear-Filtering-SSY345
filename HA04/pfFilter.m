function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, N, bResample, plotFunc)
    % PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
    % state-space model.
    %
    % Input:pfFilter(x_0
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   Y           [m x K] Measurement sequence to be filtered
    %   proc_f      Handle for process function f(x_k-1)
    %   proc_Q      [n x n] process noise covariance
    %   meas_h      Handle for measurement model function h(x_k)
    %   meas_R      [m x m] measurement noise covariance
    %   N           Number of particles
    %   bResample   boolean false - no resampling, true - resampling
    %   plotFunc    Handle for plot function that is called when a filter
    %               recursion has finished.
    % Output:
    %   xfp         [n x K] Posterior means of particle filter
    %   Pfp         [n x n x K] Posterior error covariances of particle filter
    %   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
    %   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

    
    n = size(x_0,1);
    K = size(Y,2);

    % allocate memory
    xfp = zeros(n,K);
    Pfp = zeros(n,n,K);
    Xp = zeros(n,N,K);
    Wp = zeros(N,K);
    
    % sample initial particles around prior distribution
    Xp(:,:,1) = mvnrnd(x_0,P_0,N)';
    Wp(:,1)  = 1/N * ones(1,N);
    
    for k=2:K+1
        % perform a particle filter step for the next measurement
        [Xp(:,:,k), Wp(:,k)] = pfFilterStep(Xp(:,:,k-1), Wp(:,k-1)', Y(:,k-1), proc_f, proc_Q, meas_h, meas_R);
        % plotFunc(k, X_k, W_k, xf, Pf, bResample, sigma, ax);
        % resample
        if bResample
            [Xp(:,:,k), Wp(:,k), ~] = resampl(Xp(:,:,k), Wp(:,k)'); 
        end
        % estimate mean and covariance given the particles
        xfp(:,k)   = sum( Xp(:,:,k).*Wp(:,k)' ,2 );
        Pfp(:,:,k) = Wp(:,k)'.*(Xp(:,:,k) - xfp(:,k))*(Xp(:,:,k) - xfp(:,k))';
    end
    
    % remove prior from vector
    xfp = xfp(:,2:end);
    Pfp = Pfp(:,:,2:end);
    Xp  = Xp(:,:,2:end);
    Wp  = Wp(:,2:end);
end