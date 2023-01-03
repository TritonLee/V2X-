function output = ProjGra(m_d,rho,P_max)


%% projected-gradient method for solving subproblem Q1[i]

Beta = 0.1;            % estimation accuracy
gamma = 100;
sigma = sqrt((1-Beta)/2); 
sig = 1;
alpha = 1;
dis = 10;   % distance from the ES to vehicles
R = 4;
K = 6;


ITER_MAX = 20;   % max number of outer iterations (usually converge within 30 iterations)
iter = 0:ITER_MAX;
chi = 1e-4;       % target accuracy
L = 1e-4;         % constant step size chosen by Armijo rule 
c = 1;
history.dual(1)=0;
P0 = P_max/K;
P = ones(K,1)*P0;
for i = 1:length(iter)
    
    %% update theta with gradient step 
    for k = 1:K
        temp(k) = sqrt((2^R-1)*Beta*gamma*sig/(dis^(-alpha)*P(k)*sigma^4));
        I0(k) = besseli(0,temp(k));
        E0(k) = exp((Beta*gamma+sig*(2^R-1)/(P(k)*dis^(-alpha)))/(2*sigma^2));
        
        gradient(k) = -rho*(2^R-1)*sig/(2*(1-exp(-m_d))*sigma^2*(P(k))^2)*I0(k)/E0(k);   
        L = 1e-3;
        P1(k) = P(k) - L*gradient(k);  
    
    end


    % update with projected step 

    mu0 = 1/K*(sum(P)-P_max);
    for k = 1:K 
        P_temp(k) = max((P1(k)-mu0),0);
    end
    P = P_temp; % update theta0
    
    
    % compute norm value 
    Obj_D2 = norm(P)
    dual_obj = Obj_D2;
    history.dual(i+1)=dual_obj; 
    % stopping criterion
    if i>=ITER_MAX
        if abs(history.dual(i+1)-history.dual(i))<chi % line 15 of Algorithm 2
            break;
        end
    end
end
output = P;

end

%% text of figure
% plot(iter,history.dual,'b-', 'LineWidth',2);hold on;
% xlabel('Number of PG Iterations');
% ylabel('Objective Function value of problem (26)');
% grid on;


