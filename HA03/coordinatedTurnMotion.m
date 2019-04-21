function [fx, Fx] = coordinatedTurnMotion(x, T)
    % COORDINATEDTURNMOTION calculates the predicted state using a coordinated
    % turn motion model, and also calculated the motion model Jacobian
    %
    % Input:
    %   x           [5 x 1] state vector
    %   T           [1 x 1] Sampling time
    %
    % Output:
    %   fx          [5 x 1] motion model evaluated at state x
    %   Fx          [5 x 5] motion model Jacobian evaluated at state x
    %
    % NOTE: the motion model assumes that the state vector x consist of the
    % following states:
    %   px          X-position
    %   py          Y-position
    %   v           velocity
    %   phi         heading
    %   omega       turn-rate

    fx = [x(1) + T*x(3)*cos(x(4));
          x(2) + T*x(3)*sin(x(4));
          x(3);
          x(4) + T*x(5);
          x(5)];

    % Check if the Jacobian is requested by the calling function
    if nargout > 1
        
        Fx = [1 0 T*cos(x(4)) -T*x(3)*sin(x(4)) 0;
              0 1 T*sin(x(4))  T*x(3)*cos(x(4)) 0;
              0 0 1            0                0;
              0 0 0            1                T;
              0 0 0            0                1];
    end
end