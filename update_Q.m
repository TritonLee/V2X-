%% V2X communication mode in AV lane-change scenario
function [auxi,f_value] = update_Q(Beta,gamma)

% clc;
% clear;
% close all;
% m_d = 2;
% Beta = 0.3;            % estimation accuracy
% gamma = 20;
   
sigma = sqrt((1-Beta)/2); 

f = @(be)(1-marcumq(sqrt(Beta*gamma)/sigma,sqrt(be)/sigma));
B_left = 0;                          
B_right = 100;                             
[x,fval] = goldmax(f,B_left,B_right,1e-4);
auxi = x;

f_value = fval;

end