function [solution_md, f_value] = update_D(N,T,Nx,Nu,Tsim,rho)

%clc;
%clear;
%close all;
%tic;

P_max = 100;
Beta = 0.1;            % estimation accuracy
gamma = 100;
sigma = sqrt((1-Beta)/2); 
sig = 1;
alpha = 1;
dis = 10;   % distance from the ES to vehicles
R = 4;

P_ini = P_max;
temp = (2^R-1)*sig^2/(P_ini*dis^(-alpha));
f = 1-marcumq(sqrt(Beta*gamma)/sigma,sqrt(temp)/sigma);
[marcQ] = f;
   

%% road geometry
noLane = 2; 
laneWidth = 3.7; % United States road width standard
road_right = 0;
road_center = road_right + laneWidth;
road_left = road_right + noLane * laneWidth;

d_min = 0.5;
%% Reference of trajectory for the ego vehicle
Xout = zeros(N,3);   % three dimension (including x, y and theta )
Tout = zeros(N,1); 

%% generate a straight line for reference  

for k = 1:N 
    Xout(k,1) = k*T;      % ground-truth location 
    Xout(k,2) = laneWidth+laneWidth/2;        % the objective location (trajectory of Y-coordinate) 
    Xout(k,3) = 0;
    Tout(k,1) = (k-1)*T;
end
[Nr,Nc] = size(Xout);     % Nr is the number of rows of Xout


%% Tracking a constant reference trajectory
% Initialization the vehicles
for i = 1:Nr	  
    x_l(i) = 4 + 0.5*i;   % location of leader vehicle on EL with velocity 0.5 m/s
    x_t(i) = 10 + i;  % location of target vehicle on TL with velocity 1 m/s
    x_f(i) = 0 + 0.8*i;   % location of follow vehicle on TL with velocity 0.8 m/s
end

X0 = [2 laneWidth/2 0];         % initial state: first coloum x, second coloum y, third coloum theta 

% control variable initialization
vd1 = 2;                     % longitutial leader vehicle initialize velocity 2 m/s
vd2 = 0;                     % yaw angle
x_real = zeros(Nr,Nc);       % state-space vector
x_revised = zeros(Nr,Nc);     % revised-state (considered safe margin)    
x_piao = zeros(Nr,Nc);
x_safe = zeros(Nr,1);
u_real = zeros(Nr,2);
u_piao = zeros(Nr,2);
x_real(1,:) = X0;                        % initial state x0
x_piao(1,:) = x_real(1,:) - Xout(1,:);   % objective function 
q = [1 0 0;0 0.21 0;0 0 1];              % weighted positive definite matrix
Q_cell = cell(Tsim,Tsim);
for i = 1:Tsim
    for j = 1:Tsim
        if i==j
            Q_cell{i,j}=q;
        else 
            Q_cell{i,j}=zeros(Nx,Nx);
        end 
    end
end
Q = cell2mat(Q_cell);  % use cell2mat to integrate a vector merge into a function
R = 0.1*eye(Nu*Tsim,Nu*Tsim);
W = 2*diag(ones(2,1));
Constraints = [];
%x_p = zeros(Nr,1);

%% state-space model in MPC
for i = 1:Nr
    t_d = Xout(i,3);   % yaw angle 
    a = [1    0   -vd1*sin(t_d)*T;
         0    1   vd1*cos(t_d)*T;
         0    0   1;];
    b = [cos(t_d)*T   0;
         sin(t_d)*T   0;
         0            T;];     
    A_cell = cell(Tsim,1);
    B_cell = cell(Tsim,Tsim);
     for j=1:1:Tsim
        A_cell{j,1} = a^j;
        for k = 1:1:Tsim
           if k <= j
                B_cell{j,k} = (a^(j-k))*b;  
           else
                B_cell{j,k} = zeros(Nx,Nu); 
           end
        end
    end
    A = cell2mat(A_cell);
    B = cell2mat(B_cell);
    H = 2*(B'*Q*B+R);
    H = (H+H')/2; % generate symmetric H
    f = 2*B'*Q*A*x_piao(i,:)';
    A_cons = [];
    b_cons = [];
    lb = -2*ones(Nu*Tsim,1);    % lower bound of velocity 
    ub = 2*ones(Nu*Tsim,1);     % upper bound of yaw rate 
    
        yalmip('clear')

    % Define variables
    Z = sdpvar(2,1);
    m_d = sdpvar(1,1);
    C = sdpvar(Nu*Tsim,1);
       
    Constraints = [m_d >=0,lb <= C <= ub, Z(2)==x_real(i,2)];
    if x_real(i,2) > laneWidth      % condition on moving to target lane
        Constraints = [Constraints,Z(1) + m_d <= (x_t(i) - d_min)];
        Constraints = [Constraints,Z(1) - m_d >= (x_f(i) + d_min)];
    else    % condition on moving on ego lane
        Constraints = [Constraints,Z(1) + m_d <= (x_l(i) - d_min)];
    end

    Objective = 1/2*C'*H*C + f'*C+1/2*Z'*W*Z - rho/(1-exp(-m_d))*marcQ;
    sol = optimize(Constraints,Objective);
    


    %[C,fval(i,1)] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub);  % Quadratic programming
    
    u_piao(i,1) = C(1,1);  % optimize control variable  (velocity)
    u_piao(i,2) = C(2,1);  % optimize control variable  (yaw rate)
    
    vd11 = vd1 + u_piao(i,1);  % velocity change
    vd22 = vd2 + u_piao(i,2);  % yaw angle change 
    
    ve_front(i) = 1;    % velocity of the front car
    X_front(i) = 5+ve_front(i)*i;  % location of the front car
    Y_front(i) = 3;
    
    X00 = x_real(i,:);
    % interference the dynamic of the next time slot 
    XOUT = dsolve('Dx-vd11*cos(z)=0','Dy-vd11*sin(z)=0','Dz-vd22=0','x(0)=X00(1)','y(0)=X00(2)','z(0)=X00(3)');
    t=T;  % each time slot 
    x_real(i+1,1) = eval(XOUT.x);
    x_real(i+1,2) = eval(XOUT.y);
    x_real(i+1,3) = eval(XOUT.z);
    
    %Constraints = [m_d >= 0,Z(1) + m_d <= (x_t - d_min),Z(1) - m_d >= (x_f + d_min), Z(1) + m_d <= (x_l - d_min)];    
    solution_z = value(Z);
    solution_md = value(m_d);
    solution = value(Objective);   
% 
    
    if(i<Nr)
        x_piao(i+1,:) = x_real(i+1,:) - Xout(i+1,:);
    end
    u_real(i,1) = vd1 + u_piao(i,1);
    u_real(i,2) = vd2 + u_piao(i,2);
    
    % calculate the distance for two vehicles
    d(1:Nr) = sqrt((x_real(1:Nr,1)-X_front(i)).^2 + (x_real(1:Nr,2)-Y_front(i)).^2);   

end


f_value = solution;

end