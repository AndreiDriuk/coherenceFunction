function [PhiTop, dPhiTop] = fhitop(nu_t,alpha,ztop)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
PhiTop = zeros(length(nu_t),length(ztop));
dPhiTop = zeros(length(nu_t),length(ztop));
for i=1:(length(nu_t))
  PhiTop(i,:) = hypergeom([-nu_t(i),nu_t(i)+1],1,((1+tanh(alpha*ztop))/2));
  dPhiTop(i,:) = -nu_t(i).*(nu_t(i)+1)*alpha/2.*hypergeom([-nu_t(i)+1,nu_t(i)+2],2,((1+tanh(alpha*ztop))/2))./(cosh(alpha*ztop).^2);
end
end
