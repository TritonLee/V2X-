%function [solution_md, f_value] = update_D(N,T,Nx,Nu,Tsim,rho)

clc;
clear;
%close all;
tic;

N = 10;     
T = 1;                % time slot / time interval 
Nx = 3;               % state variables 
Nu = 2;               % control variables
Tsim = 10;            % simulation time



%% road geometry
noLane = 2; 
laneWidth = 3.7; % United States road width standard
road_right = 0;
road_center = road_right + laneWidth;
road_left = road_right + noLane * laneWidth;

d_min = 5.2;
Len_car = 4.7;
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
% Initialization the surrounding vehicles

    

for i = 1:Nr	  
    x_l(i) = 30 + 1.38*i;   % location of leader vehicle on EL with velocity 5 km/h
    x_t(i) = 60 + 11.1*i;     % location of target vehicle on TL with velocity 40 km/h
    x_f(i) = 10 + 2.78*i;   % location of follow vehicle on TL with velocity 10 km/h
end
Num = 100;
Times = zeros(1,Num);
X0 = [15 laneWidth/2 0];         % initial state: first coloum x, second coloum y, third coloum theta 
pro = [0.1:0.1:0.9];            % outage probability due to imperfect CSI
    
for p = 1:length(pro)
    count = 0;
    for t = 1:length(Times)
        %I_prb(p) = pro(p)*unifrnd(-1,1);      % perturbation of parameters from -1 to 1
        I_prb(p) = -rand(1)*pro(p);
        % control variable initialization
        vd1 = 3;                     % longitutial leader vehicle initialize velocity 1 m/s
        vd2 = 0;                     % yaw angle
        x_real = zeros(Nr,Nc);       % state-space vector
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
            lb = -6*ones(Nu*Tsim,1);    % lower bound of velocity 
            ub = 6*ones(Nu*Tsim,1);     % upper bound of yaw rate 

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
            
            x_ll(i) = 30*(1+I_prb(p)) + 1.38*i;   % location of leader vehicle on EL with velocity 5 km/h
            y_ll(i) = (laneWidth/2);
            x_tt(i) = 60*(1+I_prb(p)) + 11.1*i;  % location of target vehicle on TL with velocity 40 km/h
            y_tt(i) = (laneWidth+laneWidth/2);
            x_ff(i) = 10*(1+I_prb(p)) + 2.78*i;   % location of follow vehicle on TL with velocity 10 km/h
            y_ff(i) = (laneWidth+laneWidth/2);

            d1(i) = sqrt((x_real(i,1)-x_ll(i)).^2 + (x_real(i,2)-y_ll(i)).^2);   
            d2(i) = sqrt((x_real(i,1)-x_tt(i)).^2 + (x_real(i,2)-y_tt(i)).^2);   
            d3(i) = sqrt((x_real(i,1)-x_ff(i)).^2 + (x_real(i,2)-y_ff(i)).^2);
            %count = count;
            if (d1(i) >= Len_car) && (d2(i) >= Len_car) && (d3(i) >= Len_car)
                count = count;
            else
                count = count+1;
                break
            end
            
        end
        
    end
    ratio1(p) = count/Num;   % collision probability 
    
end


toc;
save('basedline.mat','ratio1');



plot(pro,ratio1,'k-','LineWidth',2);hold on;
xlabel('Outage probability');
ylabel('Collision probability');

grid on; 
%axis([0 20 0 4.5]);



%end