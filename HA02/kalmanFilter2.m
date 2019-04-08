function [x_k_k, P_k_k, x_k_km1, P_k_km1, v, S, K] = kalmanFilter2(Y, x_0, P_0, A, Q, H, R)
    % KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
    %
    % Input:
    %   Y           [m x N] Measurement sequence
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   A           [n x n] State transition matrix
    %   Q           [n x n] Process noise covariance
    %   H           [m x n] Measurement model matrix
    %   R           [m x m] Measurement noise covariance
    %
    % Output:
    %   x           [n x N] Estimated state vector sequence
    %   P           [n x n x N] Filter error convariance
    %

    % Parameters
    N = size(Y,2);

    n = length(x_0);
    m = size(Y,1);

    % Data allocation
    x_k_k = zeros(n,   N +1);
    P_k_k = zeros(n,n, N +1);

    % filter
    x_k_k(:,1)   = x_0;
    P_k_k(:,:,1) = P_0;
    
    for i=1:N
        %% prediction step: compute p(x_k | y_1:k-1) from p(x_k-1 | y_1:k-1)
        
        x_k_km1(:,i)   = A * x_k_k(:,i);
        P_k_km1(:,:,i) = A * P_k_k(:,:,i) * A' + Q;
        
        %% update step: compute p(x_k | y_1:k) from p(x_k | y_1:k-1)
        
        % inovation mean
        v(:,i) = Y(:,i) - H * x_k_km1(:,i);
        % inovation covariance
        S(:,:,i) = H * P_k_km1(:,:,i) * H' + R;
        % kalman gain
        K(:,:,i) = P_k_km1(:,:,i) * H' / S(:,:,i);

        % updated state mean
        x_k_k(:,i+1) = x_k_km1(:,i) + K(:,:,i) * v(:,i);
        % updated state covariance
        P_k_k(:,:,i+1) = P_k_km1(:,:,i) - K(:,:,i) * S(:,:,i) * K(:,:,i)';
    end
    
    % exclude prior x0~N(x_0, P_0) from posterior
    x_k_k = x_k_k(:,2:end);
    P_k_k = P_k_k(:,:,2:end);
    
end
