function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, f, Q, h, R, N, bResample, plotFunc, pmap)
    % PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
    % state-space model.
    %
    % Input:pfFilter(x_0
    %   x_0         [n x 1] Prior mean or 
    %            OR [n x N] initial particle positions
    %   P_0         [n x n] Prior covariance
    %   Y           [m x K] Measurement sequence to be filtered
    %   f           Handle for process function f(x_k-1)
    %   Q           [n x n] process noise covariance
    %   h           Handle for measurement model function h(x_k)
    %   R           [m x m] measurement noise covariance
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
    % if x_0 has only one column, x_0 is the mean of the prior
    % otherwise x_0 are the initial particles states
    if size(x_0,2) == 1
        Xp(:,:,1) = mvnrnd(x_0,P_0,N)';
    else
        Xp(:,:,1) = x_0;
    end
    
    Wp(:,1)  = 1/N * ones(1,N);
    
    j = 1:N;
    for k=2:K+1
        
        Xp_km1 = Xp(:,:,k-1);
        Wp_km1 = Wp(:,k-1)';
        % resample
        if bResample
            [Xp_km1, Wp_km1, j] = resampl(Xp_km1, Wp_km1);        
        end
            
        % perform a particle filter step for the next measurement
        [Xp(:,:,k), Wp(:,k)] = pfFilterStep( Xp_km1, Wp_km1, Y(:,k-1), f, Q, h, R);
        % plot particles using function handle
        if ~isempty(plotFunc)
            plotFunc(k-1, Xp(:,:,k), Xp(:,:,k-1), Wp(:,k)', j);
        end
        % if handle bmap is given: update p(y_k|x(i)_k,M)=p(y_k|x(i)_k)*bmap(M|x(i)_k)
        if ~isempty(pmap)
            p_map_x = pmap(Xp(1,:,k),Xp(2,:,k));
            Wp(:,k) = Wp(:,k) .* p_map_x;
            Wp(:,k) = Wp(:,k)/sum(Wp(:,k));
        end
% % % %         % resample
% % % %         if bResample
% % % %             [Xp(:,:,k), Wp(:,k), j] = resampl(Xp(:,:,k), Wp(:,k)');
% % % %         end
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