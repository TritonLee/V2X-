%% trajectory of the baseline policy


clc;
clear;
%close all;
tic;
   
T = 1;                % time slot / time interval 
Nx = 3;               % state variables 
Nu = 2;               % control variables
Tsim = 7;            % simulation time
K = Tsim;

%% road geometry
noLane = 2; 
laneWidth = 3.72; % United States road width standard
road_right = 0;
road_center = road_right + laneWidth;
road_left = road_right + noLane * laneWidth;
d_min = 5.2;  % the length of vehicle + safe distance (0.5 m)
Len_car = 4.7;
%% Reference of trajectory for the ego vehicle
Xout = zeros(K,3);   % three dimension (including x, y and theta )
Tout = zeros(K,1); 

%% generate a straight line for reference  

for k = 1:K 
    Xout(k,1) = 23 + 1*k*T;      % ground-truth location for (X - coordinate) 
    Xout(k,2) = laneWidth+laneWidth/2; % the objective location (trajectory of Y-coordinate) 
    Xout(k,3) = 0;    %  the reference for angle raw rate
    Tout(k,1) = (k-1)*T;
end
[Nr,Nc] = size(Xout);     % Nr is the number of rows of Xout


%% Tracking a constant reference trajectory
% Initialization the surrounding vehicles

for i = 1:Nr	  
    x_l(i) = 30 + 1.38*i;   % location of leader vehicle on EL with velocity 5 km/h
    x_t(i) = 40 + 7*i;     % location of target vehicle on TL with velocity 40 km/h
    x_f(i) = 15 + 2.2*i;   % location of follow vehicle on TL with velocity 10 km/h
end

X0 = [23 laneWidth/2 0];         % initial state: first coloum x, second coloum y, third coloum theta 

%I_prb(p) = pro(p)*unifrnd(-1,1);      % perturbation of parameters from -1 to 1
I_prb = 0.9*0.1;
% control variable initialization
vd1 = 3;              % longitutial leader vehicle initialize velocity 3 m/s
vd2 = 0;              % yaw angle
x_real = zeros(Nr,Nc);      % state-space vector
x_revised = zeros(Nr,Nc);    % revised-state (considered safe margin)    
x_piao = zeros(Nr,Nc);
x_safe = zeros(Nr,1);
u_real = zeros(Nr,2);
u_piao = zeros(Nr,2);
x_real(1,:) = X0;                        % initial state x0
x_piao(1,:) = x_real(1,:) - Xout(1,:);   % objective function 
q = [1 0 0;0 0.21 0;0 0 1];               % weighted positive definite matrix
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
x_ll = zeros(Nr,1);y_ll = zeros(Nr,1);
x_tt = zeros(Nr,1);y_tt = zeros(Nr,1);
x_ff = zeros(Nr,1);y_ff = zeros(Nr,1);
d1 = zeros(Nr,1); d2 = zeros(Nr,1); d3 = zeros(Nr,1);
d_11 = zeros(Nr,1); d22 = zeros(Nr,1); d33 = zeros(Nr,1);
count = 0;
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
    lb = -6*ones(Nu*Tsim,1);    % lower bound of velocity 6m/s
    ub = 6*ones(Nu*Tsim,1);     % upper bound of velocity 

    %% use YALMAP solve convex problem (quadratic programming)
    yalmip('clear')
    % Define variables
    Z = sdpvar(2,1);   % postion of vehicle in networks (X,Y)
    C = sdpvar(Nu*Tsim,1);     % control variables (velocity and yaw rate)
    Constraints = [lb <= C <= ub];
    if x_real(i,2) > laneWidth      % condition on moving to target lane
        Constraints = [Constraints,Z(1)<= (x_t(i) - d_min)];
        Constraints = [Constraints,Z(1)>= (x_f(i) + d_min)];
    else    % condition on moving on ego lane
        Constraints = [Constraints,Z(1)<= (x_l(i) - d_min)];
    end
    Z(2) = x_real(i,2);
    Objective = 1/2*Z'*W*Z;

    %Objective = 1/2*C'*H*C + f'*C+1/2*Z'*W*Z;

    sol = optimize(Constraints,Objective);
    %[C,fval(i,1),exitflag(i,1),output(i,1)] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub);  

    [C,fval(i,1)] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub);  % Quadratic programming

    u_piao(i,1) = C(1,1);  % optimize control variable  (velocity)
    u_piao(i,2) = C(2,1);  % optimize control variable  (yaw rate)

    vd11 = vd1 + u_piao(i,1);  % velocity change
    vd22 = vd2 + u_piao(i,2);  % yaw angle change 
    X00 = x_real(i,:);
    % interference the dynamic of the next time slot 
    XOUT = dsolve('Dx-vd11*cos(z)=0','Dy-vd11*sin(z)=0','Dz-vd22=0','x(0)=X00(1)','y(0)=X00(2)','z(0)=X00(3)');
    t=T;  % each time slot 
    x_real(i+1,1) = eval(XOUT.x);
    x_real(i+1,2) = eval(XOUT.y);
    x_real(i+1,3) = eval(XOUT.z);

    if(i<Nr)
        x_piao(i+1,:) = x_real(i+1,:) - Xout(i+1,:);
    end
    u_real(i,1) = vd1 + u_piao(i,1);
    u_real(i,2) = vd2 + u_piao(i,2);

    % calculate the distance for two vehicles
    % imperfect position of surrounding vehicles 

    x_ll(i) = 30*(1+I_prb) + 1.38*i;   % location of leader vehicle on EL with velocity 5 km/h
    y_ll(i) = (laneWidth/2);
    x_tt(i) = 40*(1+I_prb) + 7*i;  % location of target vehicle on TL with velocity 40 km/h
    y_tt(i) = (laneWidth+laneWidth/2);
    x_ff(i) = 15*(1+I_prb) + 2.2*i;   % location of follow vehicle on TL with velocity 10 km/h
    y_ff(i) = (laneWidth+laneWidth/2);

    d1(i) = sqrt((x_real(i+1,1)-x_ll(i)).^2 + (x_real(i+1,2)-y_ll(i)).^2);   
    d_11(i) = abs(x_real(i+1,1)-x_ll(i));
    d2(i) = sqrt((x_real(i+1,1)-x_tt(i)).^2 + (x_real(i+1,2)-y_tt(i)).^2);   
    d22(i) = abs(x_real(i+1,1)-x_tt(i));
    d3(i) = sqrt((x_real(i+1,1)-x_ff(i)).^2 + (x_real(i+1,2)-y_ff(i)).^2);  
    d33(i) = abs(x_real(i+1,1)-x_ff(i));
    if (d1(i) >= Len_car) && (d_11(i)>= Len_car) && (d2(i) >= Len_car) && (d22(i) >= Len_car) && (d3(i) >= Len_car) && (d33(i) >= Len_car)
        count = count;
    else
        count = count+1;
        break
    end
end
toc;
%plot(Xout(1:Nr,1),Xout(1:Nr,2),'b*');hold on;
%plot(x_ll,y_ll,'ks');hold on;
%plot(x_tt,y_tt,'gs');hold on;
%plot(x_ff,y_ff,'ms');hold on;
plot(x_real(:,1),x_real(:,2),'k+');hold on;

legend({'lead vehicle trajectory','target vehicle trajectory','follow vehicle trajectory','ego vehicle trajectory of proposed policy','ego vehicle trajectory of baseline policy',},'Interpreter','latex');
xlabel('Location in X-coordinate (m)');
ylabel('Location in Y-coordinate (m)');
grid on; 

