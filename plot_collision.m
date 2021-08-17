

clc
clear;
close all;
tic;
%%  massive mimo downlink system under AM framework 
%%  minmax method to settle throughput value of Rs
load('outage_collision.mat')
load('basedline.mat','ratio1');


X = [0.7 0.5 0.3];  % x_label 
data = [ratio(4) ratio1(4); ratio(6) ratio1(6); ratio(8) ratio1(8)]; 
%Y = [ratio(3) ratio(5) ratio(7) ratio1(3) ratio1(5) ratio1(7)];
b = bar(data);
ch = get(b,'children');
%set(ch,'FaceVertexCData',[1 0 1; 0 0 0;]);
set(gca,'XTickLabel',{'0.3','0.5','0.7'});
legend('Proposed Policy','Baseline Policy');
xlabel('Outage probability');
ylabel('Collision ratio');
grid on; 
toc;