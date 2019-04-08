function [x, P] = linearPrediction(x, P, A, Q)
    % LINEARPREDICTION calculates mean and covariance of predicted state
    %   density using a linear Gaussian model.
    %
    % Input:
    %   x           [n x 1] Prior mean
    %   P           [n x n] Prior covariance
    %   A           [n x n] State transition matrix
    %   Q           [n x n] Process noise covariance
    %
    % Output:
    %   x           [n x 1] predicted state mean
    %   P           [n x n] predicted state covariance
    %

    % prediction step: compute p(x_k | y_1:k-1) from p(x_k-1 | y_1:k-1)
    %
    % Use motion model for prediction 
    % \hat{x}_{k|k-1} = A_{k-1} * \hat{x}_{k-1|k-1}
    x = A * x;
    % Compute covariance of new prediction:
    % Cov(x_{k|k-1}) = Cov(A_{k-1} * x_{k-1|k-1} + q) = A*P*A' + Q
    P = A * P * A' + Q;

end