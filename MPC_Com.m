
clc;
clear;
%close all;


%% Reference of trajectory for the ego vehicle
laneWidth = 3.72; % United States road width standard
N = 6;     
T = 1;            % time slot / time interval 
Xout = zeros(N,3);   % three dimension (including x, y and theta )
Tout = zeros(N,1); 

%% generate a straight line for reference  

for k = 1:N 
    Xout(k,1) = 12 + 3*k*T;      % ground-truth location 
    Xout(k,2) = laneWidth+laneWidth/2;           % the objective location (trajectory of Y-coordinate) 
    Xout(k,3) = 0;
    Tout(k,1) = (k-1)*T;
end
[Nr,Nc] = size(Xout);     % Nr is the number of rows of Xout 


%% Tracking a constant reference trajectory
% Initialization the vehicles

Nx = 3;               % state variables 
Nu = 2;               % control variables
Tsim = 10;            % simulation time
X0 = [15 laneWidth/2 0];          % initial state: first coloum x, second coloum y, third coloum theta 
pro = 0.3;            % outage probability due to imperfect CSI

% Mobile Robot variable Model
% control variable initialization
vd1 = 5;                     % longitutial leader vehicle initialize velocity 2 m/s
vd2 = 0;                     % yaw angle
x_real = zeros(Nr,Nc);       % state-space vector
x_piao = zeros(Nr,Nc);
u_real = zeros(Nr,2);
u_piao = zeros(Nr,2);
x_real(1,:) = X0;            % initial state x0
x_piao(1,:) = x_real(1,:) - (Xout(1,:));   % objective function 
X_PIAO = zeros(Nr,Nx*Tsim);  
q = [1 0 0;0 0.21 0;0 0 1];  % weighted matrix
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


%% state-space model in MPC
for i = 1:Nr
    t_d = Xout(i,3);
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
                B_cell{j,k} = (a^(j-k))*b;  %% product of two matrices
           else
                B_cell{j,k} = zeros(Nx,Nu);  %% linear programming
           end
        end
    end
    A = cell2mat(A_cell);
    B = cell2mat(B_cell);
    
    H = 2*(B'*Q*B+R);
    
    x_piao(i,:) = pro.*x_piao(i,:);
    f = 2*B'*Q*A*x_piao(i,:)';
    A_cons = [];
    b_cons = [];
    lb = [-1;-1];   
    ub = [1;1];
    [X,fval(i,1),exitflag(i,1),output(i,1)] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub);  
    % Quadratic programming, H positive definite matrix, f linear vector, A_cons: inequality matrix, b_cons: right vector, lb: lower bound of variable, ub: upper bound;



    u_piao(i,1) = X(1,1);
    u_piao(i,2) = X(2,1);
    %Tvec = [0:0.05:4];
    X00 = x_real(i,:);
    vd11 = vd1 + u_piao(i,1);
    vd22 = vd2 + u_piao(i,2);
    
    ve_front(i) = 0.5;    % velocity of the front car
    
    X_front(i) = 5+ve_front(i)*i;  % location of the front car
    Y_front(i) = 3;
    
    % interference the dynamic of the next time slot 
    XOUT = dsolve('Dx-vd11*cos(z)=0','Dy-vd11*sin(z)=0','Dz-vd22=0','x(0)=X00(1)','y(0)=X00(2)','z(0)=X00(3)');
    t=T; 
    x_real(i+1,1) = eval(XOUT.x);  % x-coordinate
    x_real(i+1,2) = eval(XOUT.y);  % y-coordinate
    x_real(i+1,3) = eval(XOUT.z);  % theta 
    if(i<Nr)
        x_piao(i+1,:) = x_real(i+1,:) - (Xout(i+1,:));
    end
    u_real(i,1) = vd1 + u_piao(i,1);
    u_real(i,2) = vd2 + u_piao(i,2);
    
    % calculate the distance for two vehicles
    
    %d(1:Nr) = sqrt((x_real(1:Nr,1)-X_front(i)).^2 + (x_real(1:Nr,2)-Y_front(i)).^2);
    
%     figure(1);
%     plot(Xout(1:Nr,1),Xout(1:Nr,2),'r*');hold on;
%     plot(x_real(:,1),x_real(:,2),'k*');hold on;
%     title('Trajectory Comparison');
%     legend({'Reference trajectory','Real trajectory with $p^{B}_\mathrm{out}=0.3$','Real trajectory with $p^{B}_\mathrm{out}=0.8$'},'Interpreter','latex');  ;
%     xlabel('Location in X-coordinate (m)');
%     ylabel('Location in Y-coordinate (m)');
%     grid on; 
%     figure(1)
%     plot(X_front,Y_front,'r*');hold on;
%     plot(x_real(:,1),x_real(:,2),'b^');hold on;
%     legend({'Front vehicle','Ego vehicle'},'Interpreter','latex');
%     xlabel('Location in X-coordinate (m)');
%     ylabel('Location in Y-coordinate (m)');
%     grid on; 
    
end

%% simulation part


figure()
%plot(Xout(1:Nr,1),Xout(1:Nr,2),'r*');hold on;
plot(x_real(:,1),x_real(:,2),'bs');hold on;
legend({'Reference trajectory','Real trajectory with $p^{B}_\mathrm{out}=0.3$','Real trajectory with $p^{B}_\mathrm{out}=0.8$'},'Interpreter','latex');
xlabel('Location in X-coordinate (m)');
ylabel('Location in Y-coordinate (m)');
grid on; 
%axis([0 20 0 4.5]);


% figure()
% subplot(2,1,1);
% plot(Tout(1:Nr),ve_front,'k--');hold on;   % velocity of the front car
% plot(Tout(1:Nr),(u_real(1:Nr,1)),'k');hold on; 
% legend('Front vehicle','Ego vehicle');
% %axis([0 20 0.8 2]);
% xlabel('time (s)');
% ylabel('longitudinal velocity (m/s)');
% grid on;  
% 
% subplot(2,1,2);
% plot(Tout(1:Nr),Xout(1:Nr,3),'k-.');hold on;
% plot(Tout(1:Nr),x_real(1:Nr,3),'k');hold on;
% legend('Front vehicle','Ego vehicle');hold on;
% xlabel('time (s)');
% ylabel('yaw angle \theta');
% grid on; 
% 
% 
% figure()
% plot(Tout(1:Nr),u_real(1:Nr,2),'b');hold on; 
% title('Yaw acceleration of the ego car');
% xlabel('time (s)');
% ylabel('\omega (m/s^2)');


% figure()
% plot(Tout(1:Nr),d(1:Nr),'r-','LineWidth',2);hold on;
% xlabel('sampling time (s)');
% ylabel('Distance between two vehicles (m)');
% grid on;

