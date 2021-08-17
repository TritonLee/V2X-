function [ value ] = Threshold (x,alpha,sigma,delta,epsilon)
% Bisect get threhold value of throughput 
% Based on S-ADP transmission scheme (fixed communication rate Rb)
% Detailed explanation goes here
    left = marcumq(sqrt(alpha)/sigma,sqrt(x)/sigma);
    value = left - (1-delta);
end
