

clc
clear;
close all;
tic;
%%  massive mimo downlink system under AM framework 
%%  minmax method to settle throughput value of Rs
load('basedline.mat','ratio1');

ratio = [0 0 0 0 0.05 0.05 0.06 0.1];
ratio2 = [0.01 0.1 0.1 0.1 0.2 0.15 0.2 0.2];
ratio1 = [0.1 0.6 0.6 0.8 0.8 1 1 1]



X = [10 20 40];  % x_label 
data = [ratio(1) ratio2(1) ratio1(1); ratio(2) ratio2(2) ratio1(2); ratio(6) ratio2(6) ratio1(6)]; 
%Y = [ratio(3) ratio(5) ratio(7) ratio1(3) ratio1(5) ratio1(7)];
b = bar(data);
ch = get(b,'children');
%set(ch,'FaceVertexCData',[1 0 1; 0 0 0;]);
set(gca,'XTickLabel',{'10','20','40'});
legend('Proposed policy','Baseline policy with communication delay','Baseline policy ignoring communication delay');
xlabel('Car speed (km/h)');
ylabel('Collision ratio');
grid on; 
toc;