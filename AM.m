%% Outer iteratively optimize the sum rate 
% output the converged objective function value, variable Xi and Theta 

%function [xi,theta_temp,fvalue] = AM(N,T,Nx,Nu,Tsim,rho)

clc;
clear;
close;

tic;
%% initialize the updated variables
Beta = 0.5;            % estimation accuracy
gamma = 10;
N = 1;     
T = 1;            % time slot / time interval 
Nx = 3;               % state variables 
Nu = 2;               % control variables
Tsim = 1;            % simulation time
rho = 1e2;



Tol = 10^(-3);                   % search tolerance for iterative algorithm
%% iteratively optimize the based on AM framework
% update xi
[md, fvalue_l] = update_D(N,T,Nx,Nu,Tsim,rho);
% update theta

[fvalue_r] = update_Q(md,Beta,gamma);
Ra_left = fvalue_l;
Ra_right = fvalue_r;


ITER_MAX=5; % max number of outer iterations (usually converge within 5 iterations)
history.dual(1)=0;
iter = 0:ITER_MAX;
for i =1:length(iter)
    %stopping criterion
    Rate_opt = (Ra_right+Ra_left)/2;
    dual_obj = Rate_opt;
    history.dual(i+1)=dual_obj; 
    if i>=ITER_MAX
        if abs(history.dual(i+1)-history.dual(i))<Tol % line 15 of Algorithm 2
            break;
        end
    end
    [md_t, fvalue_l_t] = update_D(N,T,Nx,Nu,Tsim,rho);
    [fvalue_r_t] = update_Q(md_t,Beta,gamma);


    Ra_left = fvalue_l_t;
    Ra_right = fvalue_r_t;
end
Rate_opt = (Ra_right+Ra_left)/2;
fvalue = Rate_opt;
toc;
%end

plot(iter,history.dual,'b-', 'LineWidth',2);hold on; 
% text of figure
xlabel('Number of Outer Iterations');
ylabel('Objective function value of (26)');
legend('M = 2','M = 5','M = 10');
grid on;
