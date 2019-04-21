function Y = genLinearMeasurementSequence(X, H, R)
    % GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
    % sequence X using a linear measurement model. Measurement noise is assumed to be 
    % zero mean and Gaussian.
    %
    % Input:
    %   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
    %   H           [m x n] Measurement matrix
    %   R           [m x m] Measurement noise covariance
    %
    % Output:
    %   Y           [m x N] Measurement sequence
    %

    m = size(H,1);
    % state sequence includes x0, which does not generate an observation
    N = size(X,2) -1;
    Y = zeros(m,N);
    
    % iterate to generate N samples
    for i=1:N
        % Measurement model: Y{k} = H*X{k} + r{k},   where r{k} ~ N(0,R)
        Y(:,i) = H * X(:,i+1) + mvnrnd(zeros(m,1), R)';
    end
    
end