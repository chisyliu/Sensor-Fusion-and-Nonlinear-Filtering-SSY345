function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
    % DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
    % s1 and s2, to the position given by the state vector x. Also returns the
    % Jacobian of the model at x.
    %
    % Input:
    %   x           [n x 1] State vector, the two first element are 2D position
    %   s1          [2 x 1] Sensor position (2D) for sensor 1
    %   s2          [2 x 1] Sensor position (2D) for sensor 2
    %
    % Output:
    %   hx          [2 x 1] measurement vector
    %   Hx          [2 x n] measurement model Jacobian
    %
    % NOTE: the measurement model assumes that in the state vector x, the first
    % two states are X-position and Y-position.

    
    % Procedure to calculate hx and Hx using symbolic toolbox
    
% % %     % states
% % %     syms px py real
% % %     x = [px py].'
% % %     
% % %     % sensor positions
% % %     s1 = sym('s1_',[2;1],'real')
% % %     s2 = sym('s2_',[2;1],'real')
% % %     
% % %     % sensor readings
% % %     ang1 = atan2( py-s1(2), px-s1(1) )
% % %     ang2 = atan2( py-s2(2), px-s2(1) )
% % %     
% % %     hx = [ang1; ang2]
% % %     Hx = jacobian(hx,x)

    % initialize outputs with correct sizes
    n  = size(x,1);
    hx = zeros(2,1);
    Hx = zeros(2,n);

    % calculate readings from the two sensors
    ang1 = atan2( x(2)-s1(2), x(1)-s1(1) );
    ang2 = atan2( x(2)-s2(2), x(1)-s2(1) );
    
    % output is the concatenation of the two readings
    hx(1:2,1) = [ang1; 
                 ang2];
    
    % jacobian of hx, as calculated using the symbolic toolbox
    Hx(1:2,1:2) = [-(x(2)-s1(2)) / ( (x(1)-s1(1))^2 + (x(2)-s1(2))^2 ), (x(1)-s1(1)) / ((x(1)-s1(1))^2 + (x(2)-s1(2))^2);
                   -(x(2)-s2(2)) / ( (x(1)-s2(1))^2 + (x(2)-s2(2))^2 ), (x(1)-s2(1)) / ((x(1)-s2(1))^2 + (x(2)-s2(2))^2)];
end



