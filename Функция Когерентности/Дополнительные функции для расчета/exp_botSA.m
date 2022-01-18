function [ EXP_bot] = exp_botSA( sigma2,le,lb,kd,kplm,ks,alpha, z_bot )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%A =ones(length(k2),length(z_top));
EXP_bot=exp(-sigma2*le*kd.^2*kplm^4./(8*(ks.^2-1/4*kd.^2).^2)*(1/alpha+z_bot.*(1-z_bot.^2/(3*lb^2))));
end
