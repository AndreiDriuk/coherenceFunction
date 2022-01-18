function [ PhiBot, DPhiBot ] = phibot(DEV,DOV,DDEV,DDOV,PhiTop,dPhiTop,kd,z_bot)
[~,fk]=size(PhiTop);
[bw,bk] = size(DEV);
A=PhiTop(:,fk)./DEV(:,1);
C=dPhiTop(:,fk)./DDOV(:,1);
PhiBot = ones(bw,bk);
DPhiBot = ones(bw,bk);
for i=1:bk
PhiBot(:,i) = A.*DEV(:,i)+C.*DOV(:,i);
DPhiBot(:,i) = A.*DDEV(:,i)+C.*DDOV(:,i);
end
for j=1:bw
   if kd(j)==0 
      PhiBot(j,:) = A(j,:).*DEV(j,:)+dPhiTop(j,fk)./PhiTop(j,fk).*z_bot;
      DPhiBot(j,:) = A(j,:).*DDEV(j,:)+dPhiTop(j,fk)./PhiTop(j,fk).*DDOV(j,:);
   end
end