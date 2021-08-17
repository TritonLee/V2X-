clc
clear;
%close all;
tic;

%%  Imperfect of ADP method 
%% variable-- CSI gamma_b

T_b = 10^(0/10);                  % parameter related to sigmal_B
T_e = 10^(-5/10);                 % parameter related to sigmal_E
P_max = 10^(10/10);               % maximum SNR 
epsilon = 0.1;                    % the preset SOP value 
delta = 0.5;                      % the preset ROP value
alpha = 0.3;                      % estimation accuracy 
sigma = sqrt((1-alpha)*T_b/2); 

% calculte the threshold value miu
miu_left = 0;                     
miu_right = 50;                  
miu_temp = (miu_left + miu_right)/2;          
miu_low = Threshold (miu_left,alpha,T_e,sigma,delta,epsilon); 
miu_up = Threshold (miu_right,alpha,T_e,sigma,delta,epsilon);   
while (abs(miu_right - miu_left)>1e-5)
    miu_temp = (miu_right+miu_left)/2;         
    Miu_temp = Threshold (miu_temp,alpha,T_e,sigma,delta,epsilon);
    if (Miu_temp*miu_low>0)
        miu_left = miu_temp;
        miu_low = Miu_temp;
    elseif (Miu_temp*miu_up>0)
        miu_right = miu_temp;
        miu_up = Miu_temp;
    else
        break
    end
end
 

% the optimal threshold

if delta > 1-exp(T_e*log(epsilon)/(2*(sigma)^2))
    miu_opt = 0.01;

else
    miu_opt = miu_temp;

end

delta_gamma = 0.1;
%gamma_dB = 10*log10(miu_opt):delta_gamma:10;  % range of X axis
gamma = 10.^(10./10);
Num = length(gamma);

for ii = 1:length(gamma)
        threshold(ii) = -T_e*log(epsilon)/gamma(ii);       
        % bisection method to find the upper bound of beta parameter  
        
        beta_up_left = threshold(ii);                             
        beta_up_right = 300;                                       
        beta_up_temp = (beta_up_left + beta_up_right)/2;          
        beta_up_low = beta1(beta_up_left,alpha,sigma,gamma(ii),delta); 
        beta_up_up = beta1(beta_up_right,alpha,sigma,gamma(ii),delta);   
        while (abs(beta_up_right-beta_up_left)>1e-7)
            beta_up_temp = (beta_up_left + beta_up_right)/2;        
            Beta_up_temp = beta1(beta_up_temp,alpha,sigma,gamma(ii),delta);
            if (Beta_up_temp*beta_up_low>0)
                beta_up_left = beta_up_temp;
                beta_up_low = Beta_up_temp;
            elseif (Beta_up_temp*beta_up_up>0)
                beta_up_right = beta_up_temp;
                beta_up_up = Beta_up_temp;
            else
                break
            end
        end
        beta_up_opt(ii) = beta_up_temp;
        Be_up(ii) = beta_up_opt(ii); 
        
        % determine optimization value of beta with GSM method       
        f = @(Be)log2((1+P_max*Be)/(1-P_max*T_e*log(epsilon)))*marcumq(sqrt(alpha*gamma(ii))/sigma,sqrt(Be)/sigma);
        B_left = threshold(ii);                    
        B_right = Be_up(ii);                            
        [x,fval] = goldmax(f,B_left,B_right,1e-6);
        B_opt(ii) = x;    % optimal value of beta
        
        f_gamma(ii) = exp(-gamma(ii)/T_b)/T_b; 
        Rs(ii) = log2((1+P_max*B_opt(ii))/(1-P_max*T_e*log(epsilon)));
        Rrs(ii) = Rs(ii)*marcumq(sqrt(alpha*gamma(ii))/sigma,sqrt(B_opt(ii))/sigma);
        Rb(ii) = log2(1+P_max*B_opt(ii));    % optimal communication rate
        pro(ii) = marcumq(sqrt(alpha*gamma(ii))/sigma,sqrt(B_opt(ii))/sigma);
        
        
        ex_po(ii)= exp((alpha*gamma(ii)+B_opt(ii))/(2*sigma^2));
        Bess(ii) = besseli(0,sqrt(alpha*gamma(ii)*B_opt(ii))/sigma^2);
        
        th_ta(ii) = marcumq(sqrt(alpha*gamma(ii))/sigma,sqrt(B_opt(ii))/sigma)*ex_po(ii)/Bess(ii);
                    
        if gamma(ii)<miu_opt
            Rb(ii) = 0;
            Rs(ii) = 0;
            Rrs(ii) = 0;
        else
            Rb(ii) = Rb(ii);
            Rs(ii) = Rs(ii);
            Rrs(ii) = Rrs(ii);
        end
        
        P_opt(ii) = P_max;

end
toc;

plot(gamma_dB,B_opt,'r-','LineWidth',2);hold on;    % communication rate R_b
%plot(gamma_dB,Be_up,'b-','LineWidth',2);hold on;    % communication rate R_b
plot(gamma_dB,pro,'k-','LineWidth',2);hold on;    % communication rate R_b

%plot(gamma_dB,Rb,'bs--','LineWidth',1.5);hold on;    % communication rate R_b
%plot(gamma_dB,Rs,'r-','LineWidth',1.5);hold on;    % secrecy rate R_s
%plot(gamma_dB,Rrs,'k-','LineWidth',1.5);hold on;  % reliable secrecy rate R_rs

legend off;
%% text of figure
xlabel({'$\hat{\gamma}_\mathrm{B}$'},'Interpreter','latex');
ylabel('Rate parameters (bits/s/Hz)');


legend({'$R_\mathrm{B}, \overline{R}_\mathrm{B} = 4$','$R_\mathrm{s},\overline{R}_\mathrm{B} = 4$','$R_\mathrm{B},\overline{R}_\mathrm{B} = 3$',...
    '$R_\mathrm{s},\overline{R}_\mathrm{B} = 3$','$R_\mathrm{B},\overline{R}_\mathrm{B} = \infty$','$R_\mathrm{s},\overline{R}_\mathrm{B} = \infty$'},'Interpreter','latex');

legend off;
axis([-14 -4 0.2 0.9]);
grid on;


legend({'Proposed $N = 100$','Ignoring CSI uncertainty $N = 100$','Proposed $N = 20$','Ignoring CSI uncertainty $N = 20$'},'Interpreter','latex');
legend({'Proposed $\delta = 0.3$','Ignoring CSI uncertainty $\delta = 0.3$','Proposed $\delta = 0.5$','Ignoring CSI uncertainty $\delta = 0.5$'},'Interpreter','latex');


legend({'Proposed $P_\mathrm{max} = 0$dBm','Ignoring CSI uncertainty $P_\mathrm{max} = 0$dBm','Proposed $P_\mathrm{max} = 10$dBm','Ignoring CSI uncertainty $P_\mathrm{max} = 10$dBm'},'Interpreter','latex');

legend({'Proposed Solution','Globally Optimal Solution'},'Interpreter','latex');


legend({'Proposed solution','Ignoring CSI uncertainty','FPS scheme','RPS scheme'},'Interpreter','latex');



legend({'Proposed $M = 5$','TDMA $M = 5$','Proposed $M = 10$','TDMA $M = 10$'},'Interpreter','latex');

legend({'Proposed $J = 2$','TDMA $J = 2$','Proposed $J = 10$','TDMA $J = 10$'},'Interpreter','latex');

legend({'$P = 10$dB','$P = 0$dB','$P = -5$dB'},'Interpreter','latex');


legend({'$M = 2$','$M = 5$','$M = 10$'},'Interpreter','latex');

ylabel('Security guaranteed sum-rate (bits/s/Hz)');

ylabel('Average secure EE (bits/J)');


ylabel('Average computation time (s)');

xlabel({'$P_\mathrm{max}$ (dBm)'},'Interpreter','latex');

xlabel({'$\varepsilon_k$'},'Interpreter','latex');


xlabel({'Number of Users in Each Cluster $K_m$'},'Interpreter','latex');
xlabel({'$P_\mathrm{max}$ (dBm)'},'Interpreter','latex');
ylabel({'$g(\Xi^{(m,j)}_k,y_k)$'},'Interpreter','latex');
ylabel({'$f(\{\Theta^{(m,j)}_k\}_{k=1}^{K_m})$'},'Interpreter','latex');
ylabel({'Objective function value of $\mathcal{P}2^{[m,j]}$'},'Interpreter','latex');

%text(gamm_cri, 0.2, {'$\varsigma$'},'Interpreter','latex','FontSize', 18);
text({'$\mu$'},'Interpreter','latex','FontSize', 18);
ylabel({'$M$'},'Interpreter','latex');

xlabel({'$M$'},'Interpreter','latex');

Security guaranteed sum-rate



