function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
    % NONLINKFUPDATE calculates mean and covariance of predicted state
    %   density using a non-linear Gaussian model.
    %
    % Input:
    %   x           [n x 1] Predicted mean
    %   P           [n x n] Predicted covariance
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
    
    n = size(x,1);
    switch type
        case 'EKF'
            % calculate h(x) and jacobian(h(x),x)
            [hx, dhx] = h(x);
            % calculate inovation covariance and kalman gain
            S = dhx * P * dhx' + R;
            K = P * dhx' / S;
            % update
            x = x + K * ( y - hx );
            P = P - K * S * K';
            
        case 'UKF'
            % compute sigma points
            [SP,W] = sigmaPoints(x, P, type);
            
            % estimate desired moments
            yhat = 0*y;
            for i=1:2*n+1
                yhat = yhat + h(SP(:,i)) * W(i);
            end
            Pxy = 0*x*y';
            for i=1:2*n+1
                Pxy = Pxy + (SP(:,i)-x)*(h(SP(:,i))-yhat)' * W(i);
            end
            S = R;
            for i=1:2*n+1
                S = S + (h(SP(:,i))-yhat)*(h(SP(:,i))-yhat)' * W(i);
            end
            
            % calculate updated mean and covariance
            x = x + Pxy / S * ( y - yhat );
            P = P - Pxy / S * Pxy';
            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
            
        case 'CKF'
            % compute sigma points
            [SP,W] = sigmaPoints(x, P, type);
            
            % estimate desired moments
            yhat = 0*y;
            for i=1:2*n
                yhat = yhat + h(SP(:,i)) * W(i);
            end
            Pxy = 0*x*y';
            for i=1:2*n
                Pxy = Pxy + (SP(:,i)-x)*(h(SP(:,i))-yhat).' * W(i);
            end
            S = R;
            for i=1:2*n
                S = S + (h(SP(:,i))-yhat)*(h(SP(:,i))-yhat).' * W(i);
            end
            
            % calculate updated mean and covariance
            x = x + Pxy / S * ( y - yhat );
            P = P - Pxy / S * Pxy';
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end
end