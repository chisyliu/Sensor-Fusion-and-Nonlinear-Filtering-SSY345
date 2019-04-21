function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
    % NONLINKFUPDATE calculates mean and covariance of predicted state
    %   density using a non-linear Gaussian model.
    %
    % Input:
    %   x           [n x 1] Prior mean
    %   P           [n x n] Prior covariance
    %   y           [m x 1] measurement vector
    %   h           Measurement model function handle
    %               [hx,Hx]=h(x) 
    %               Takes as input x (state), 
    %               Returns hx and Hx, measurement model and Jacobian evaluated at x
    %               Function must include all model parameters for the particular model, 
    %               such as sensor position for some models.
    %   R           [m x m] Measurement noise covariance
    %   type        String that specifies the type of non-linear filter
    %
    % Output:
    %   x           [n x 1] updated state mean
    %   P           [n x n] updated state covariance
    %

    switch type
        case 'EKF'
            [hx, dhx] = h(x);
            
            S = dhx * P * dhx' + R;
            K = P * dhx' / S;
            
            x = x + K * ( y - hx );
            P = P - K * S * K';
            
            
        case 'UKF'
    
            % Your UKF update here
            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
            
        case 'CKF'
    
            % Your CKF update here
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

