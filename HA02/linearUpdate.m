function [x, P] = linearUpdate(x, P, y, H, R)
    % LINEARUPDATE calculates mean and covariance of predicted state
    %   density using a linear Gaussian model.
    %
    % Input:
    %   x           [n x 1] Prior mean
    %   P           [n x n] Prior covariance
    %   y           [m x 1] Measurement
    %   H           [m x n] Measurement model matrix
    %   R           [m x m] Measurement noise covariance
    %
    % Output:
    %   x           [n x 1] updated state mean
    %   P           [n x n] updated state covariance
    %

    % update step: compute p(x_k | y_1:k) from p(x_k | y_1:k-1)
    %
    % inovation mean
    v = y - H * x;
    % inovation covariance
    S = H * P * H' + R;
    % kalman gain
    K = P * H' / S;

    % updated state mean
    x = x + K * v;
    % updated state covariance
    P = P - K * S * K';

end