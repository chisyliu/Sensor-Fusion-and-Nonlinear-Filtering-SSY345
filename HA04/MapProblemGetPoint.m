
clear all
%This file draw a map of the village, and allows us to manually
%draw the trajectory of the vehicle.

figure(1)
clf
hold on
plot([1+i 1+9*i 5+9*i])
plot([7+9*i 11+9*i 11+i 7+i]);plot([5+i 1+i])
plot([2+5.2*i 2+8.3*i 4+8.3*i 4+5.2*i 2+5.2*i])%House 1
plot([2+3.7*i 2+4.4*i 4+4.4*i 4+3.7*i 2+3.7*i])%House 2
plot([2+2*i 2+3.2*i 4+3.2*i 4+2*i 2+2*i])%House 3
plot([5+i 5+2.2*i 7+2.2*i 7+i])%House 4
plot([5+2.8*i 5+5.5*i 7+5.5*i 7+2.8*i 5+2.8*i])%House 5
plot([5+6.2*i 5+9*i]);plot([7+9*i 7+6.2*i 5+6.2*i])%House 6
plot([8+4.6*i 8+8.4*i 10+8.4*i 10+4.6*i 8+4.6*i])%House 7
plot([8+2.4*i 8+4*i 10+4*i 10+2.4*i 8+2.4*i])%House 8
plot([8+1.7*i 8+1.8*i 10+1.8*i 10+1.7*i 8+1.7*i])%House 9

axis([0.8 11.2 0.8 9.2])
title('A map of the village','FontSize',20)

disp('Start clicking in the village to create a trajectory!')
disp('Press "Return" to finish.')

[X,Y]=ginput;
plot([X+Y*i],'-*')

Xk = [X';Y']
save Xk
