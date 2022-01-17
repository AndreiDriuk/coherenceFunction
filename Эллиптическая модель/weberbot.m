function [ DEV, DOV,DDEV,DDOV ] = weberbot(nu_b,ksi,dksi)

[l1,l2]=size(ksi);
B1=ones(l1,l2);
B2=ones(l1,l2);
DEV=ones(l1,l2);
DOV=ones(l1,l2);
F1 = ones(l1, l2);
dF1 = ones(l1, l2);
F2 = ones(l1, l2);
dF2 = ones(l1, l2);
parfor i=1:l1
    DEV(i,:)=exp(-ksi(i,:).^2/4).*hypergeom(-nu_b(i)/2,1/2,ksi(i,:).^2/2);
    DOV(i,:)=ksi(i,:).*exp(-ksi(i,:).^2/4).*hypergeom((1-nu_b(i))/2,3/2,ksi(i,:).^2/2);
    F1(i, :) = hypergeom((1-nu_b(i))/2,3/2,(ksi(i,:).^2)/2); 
    dF1(i, :) = hypergeom((3-nu_b(i))/2,5/2,(ksi(i,:).^2)/2);
    B1(i,:)=F1(i,:)-(ksi(i,:).^2).*F1(i,:)/2+((ksi(i,:).^2).*((1-nu_b(i))/3)).*dF1(i,:);
end
clear F1 F2
DDOV=dksi.*B1.*exp(-ksi.^2/4);

A=exp(-(ksi/2).^2);
parfor i=1:l1
  F2(i,:) = hypergeom(-nu_b(i)/2,1/2,(ksi(i,:).^2/2));
  dF2(i,:) = hypergeom((2-nu_b(i))/2,3/2,(ksi(i,:).^2)/2);
  B2(i,:)=ksi(i,:).*(-1/2*F2(i,:)-nu_b(i).*dF2(i,:));
end

DDEV=dksi.*A.*B2;
