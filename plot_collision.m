

clc
clear;
close all;
tic;
%%  massive mimo downlink system under AM framework 
%%  minmax method to settle throughput value of Rs

load('basedline.mat','ratio1');
ratio = [0 0 0 0.01 0.01 0.1 0.1 0.1];
ratio2 = [0.05 0.08 0.1 0.1 0.2 0.1 0.1 0.2];
ratio1 = [0.1 0.6 0.6 0.8 0.8 1 1 1]

X = [0.1 0.3 0.5];  % x_label 
data = [ratio(1) ratio2(1) ratio1(1); ratio(3) ratio2(3) ratio1(3); ratio(5) ratio2(5) ratio1(5)]; 
%Y = [ratio(3) ratio(5) ratio(7) ratio1(3) ratio1(5) ratio1(7)];
b = bar(data);
ch = get(b,'children');
%set(ch,'FaceVertexCData',[1 0 1; 0 0 0;]);
set(gca,'XTickLabel',{'0.1','0.3','0.5'});
legend('Proposed policy','Baseline policy with communication delay','Baseline policy ignoring communication delay');
xlabel('Outage probability');
ylabel('Collision ratio');
grid on; 
toc;