function [ EXP_top] = exp_topSA( sigma2,le,kd,kplm,ks,alpha, z_top )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%A =ones(length(k2),length(z_top));
EXP_top=exp(-sigma2*le*kd.^2*kplm^4./(8*alpha*(ks.^2-1/4*kd.^2).^2)*(1+tanh(alpha*z_top)));
end

