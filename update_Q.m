%% V2X communication mode in AV lane-change scenario
function [f_value] = update_Q(P_max,Beta,gamma)

% clc;
% clear;
% close all;
% m_d = 2;
% Beta = 0.3;            % estimation accuracy
% gamma = 20;
  
sigma = sqrt((1-Beta)/2); 
sig = 1;
alpha = 2.2;
dis = 10;   % distance from the ES to vehicles
R = 4;

P_ini = P_max;
temp = (2^R-1)*sig^2/(P_ini*dis^(-alpha));
f = 1-marcumq(sqrt(Beta*gamma)/sigma,sqrt(temp)/sigma);
f_value = f;

end