function X = genNonLinearStateSequence(x_0, P_0, f, Q, N)
    % GENNONLINEARSTATESEQUENCE generates an N+1-long sequence of states using a 
    %    Gaussian prior and a linear Gaussian process model
    %
    % Input:
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   f           Motion model function handle
    %               [fx,Fx]=f(x) 
    %               Takes as input x (state), 
    %               Returns fx and Fx, motion model and Jacobian evaluated at x
    %               All other model parameters, such as sample time T,
    %               must be included in the function
    %   Q           [n x n] Process noise covariance
    %   N           [1 x 1] Number of states to generate
    %
    % Output:
    %   X           [n x N+1] State vector sequence

    n = length(x_0);
    X = zeros(n,N);
    
    % sample initial state from the prior distribution x0~N(x_0,P_0)
    X(:,1) = mvnrnd(x_0, P_0)';
    % iterate to generate N samples
    for i=1:N
        % Motion model: X{k} = f(X{k-1}) + q{k-1},   where q{k-1} ~ N(0,Q)
        X(:,i+1) = f(X(:,i)) + mvnrnd(zeros(n,1), Q)';
    end
end