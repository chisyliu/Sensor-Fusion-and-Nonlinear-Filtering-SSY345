function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
                                     Ps_kplus1, ...
                                     xf_k, ... 
                                     Pf_k, ...
                                     xp_kplus1, ...
                                     Pp_kplus1, ...
                                     f, ...
                                     T, ...
                                     sigmaPoints, ...
                                     type)
    % NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
    % density, using a non-linear Gaussian model.
    %
    % Input:
    %   xs_kplus1   Smooting estimate for state at time k+1
    %   Ps_kplus1   Smoothing error covariance for state at time k+1
    %   xf_k        Filter estimate for state at time k
    %   Pf_k        Filter error covariance for state at time k
    %   xp_kplus1   Prediction estimate for state at time k+1
    %   Pp_kplus1   Prediction error covariance for state at time k+1
    %   f           Motion model function handle
    %   T           Sampling time
    %   sigmaPoints Handle to function that generates sigma points.
    %   type        String that specifies type of non-linear filter/smoother
    %
    % Output:
    %   xs          Smoothed estimate of state at time k
    %   Ps          Smoothed error convariance for state at time k


    %% calculate Pkkp1 = P_{k,k+1|k}
    
    if strcmp( type, 'EKF' )   
        % evaluate motion model at ^x_{k|k}
        [~, df_xk] = f(xf_k,T);
        % approximate P_{k,k+1|k}
        Pkkp1 = Pf_k * df_xk';
        
    else % UKF or CKF
        % calcualte sigma points of p(x_{k|k})
        [SP,W] = sigmaPoints(xf_k, Pf_k, type);
        % predict covariance Pkkp1 = P_{k,k+1|k}
        Pkkp1 = zeros(size(xf_k,1));
        for i=1:numel(W)
            Pkkp1 = Pkkp1 + (SP(:,i)-xf_k)*(f(SP(:,i),T)-xp_kplus1).' * W(i);
        end        
    end

    %% backward recursion

    % calculate smoothing gain
    Gk =  Pkkp1 / Pp_kplus1;
    % calculate Smoothed estimate and covariance at time k
    xs = xf_k + Gk * ( xs_kplus1 - xp_kplus1 );
    Ps = Pf_k - Gk * ( Pp_kplus1 - Ps_kplus1 ) * Gk';
    
end
