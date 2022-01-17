clear 
ztop = -1000000:10000:0;
lt = 154154;
lb = 147536;
lend = 422000;
kplm  = 0.22417;
phi = 90/180*pi;
zbot = [0:10000:lb,lb];
zunder = [lb:50000:lend,lend];
z_all = [ztop,zbot,zunder];

fs=10^9;
fd = ([-6:0.2:6]*10^8)';
sigma2=0.01;

v = 299792458; 
ks = fs*2*pi/v;
kd = fd*2*pi/v;
ax = 1000;
by = ax;
le = ax;
x = -1000:100:1000;
y = -1000:100:1000;
nu_tSAx=-1/2-1/2*sqrt(1+2*1i*kd*sigma2*kplm^4*le./((ks.^2-1/4*kd.^2).^2/lt^2*ax^2));
nu_tSAy=-1/2-1/2*sqrt(1+2*1i*kd*sigma2*kplm^4*le./((ks.^2-1/4*kd.^2).^2/lt^2*by^2));
nu_bx=sqrt(sigma2)*kplm^2*lb*sqrt(1i*kd*le/(2*ax^2))./(2*(ks.^2-1/4*kd.^2))-1/2;
nu_by=sqrt(sigma2)*kplm^2*lb*sqrt(1i*kd*le/(2*by^2))./(2*(ks.^2-1/4*kd.^2))-1/2;


ksi_bx=ones(length(kd),length(zbot));
ksi_by=ones(length(kd),length(zbot));
dksi_bx=ones(length(kd),length(zbot));
dksi_by=ones(length(kd),length(zbot));

for i=1:length(kd)
    ksi_bx(i,:)=zbot.*sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(i)*le/(2*ax^2))./((ks.^2-1/4*kd(i).^2)*lb));
    ksi_by(i,:)=zbot.*sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(i)*le/(2*by^2))./((ks.^2-1/4*kd(i).^2)*lb));
    dksi_bx(i,:)=sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(i)*le/(2*ax^2))./((ks.^2-1/4*kd(i).^2)*lb));
    dksi_by(i,:)=sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(i)*le/(2*by^2))./((ks.^2-1/4*kd(i).^2)*lb));
end



[ PsiTopSAx, dPsiTopSAx ] = fhitop(nu_tSAx,1/lt,ztop);
[ PsiTopSAy, dPsiTopSAy ] = fhitop(nu_tSAy,1/lt,ztop);
GTopFunSAx = ones(1,1,length(kd),length(ztop));
GTopFunSAy = ones(1,1,length(kd),length(ztop));
for i =1:length(kd)
    if kd(i)~=0
      GTopFunSAx(1,1,i,:) = 1i*ax^2*(ks.^2-1/4*kd(i).^2)*dPsiTopSAx(i,:)./(2*kd(i)*PsiTopSAx(i,:));
      GTopFunSAy(1,1,i,:) = 1i*by^2*(ks.^2-1/4*kd(i).^2)*dPsiTopSAy(i,:)./(2*kd(i)*PsiTopSAy(i,:));
    else
      hyper=hypergeom([-nu_tSAx(i)+1,nu_tSAx(i)+2],2,((1+tanh(ztop/lt))/2));
      GTopFunSAx(1,1,i,:) = sigma2*kplm^4*le.*hyper./(8*(ks.^2-1/4*kd(i).^2)/lt*cosh(ztop/lt).^2.*PsiTopSAx(i,:));
      hyper=hypergeom([-nu_tSAy(i)+1,nu_tSAy(i)+2],2,((1+tanh(ztop/lt))/2));
      GTopFunSAy(1,1,i,:) = sigma2*kplm^4*le.*hyper./(8*(ks.^2-1/4*kd(i).^2)/lt*cosh(ztop/lt).^2.*PsiTopSAy(i,:));      
    end
end
clear hyper


FTopFunSAx=1./sqrt(PsiTopSAx);
FTopFunSAy=1./sqrt(PsiTopSAy);
FTopFunSA(1,1,:,:)=FTopFunSAx.*FTopFunSAy;

[ DEV_SAx, DOV_SAx,DDEV_SAx,DDOV_SAx ] = weberbot(nu_bx,ksi_bx,dksi_bx);
[ DEV_SAy, DOV_SAy,DDEV_SAy,DDOV_SAy ] = weberbot(nu_by,ksi_by,dksi_by);


[ PsiBotSAx, dPsiBotSAx ] = phibot(DEV_SAx,DOV_SAx,DDEV_SAx,DDOV_SAx,PsiTopSAx,dPsiTopSAx,kd,zbot);
[ PsiBotSAy, dPsiBotSAy ] = phibot(DEV_SAy,DOV_SAy,DDEV_SAy,DDOV_SAy,PsiTopSAy,dPsiTopSAy,kd,zbot);

clear DEV_SAx DOV_SAx DDEV_SAx DDOV_SAx DEV_SAy DOV_SAy DDEV_SAy DDOV_SAy
GBotFunSAx = ones(1,1,length(kd),length(zbot));
GBotFunSAy = ones(1,1,length(kd),length(zbot));

for i =1:length(kd)
    if kd(i)~=0
      GBotFunSAx(1,1,i,:) = 1i*(ks.^2-1/4*kd(i).^2)*ax^2*dPsiBotSAx(i,:)./(2*kd(i)*PsiBotSAx(i,:));
      GBotFunSAy(1,1,i,:) = 1i*(ks.^2-1/4*kd(i).^2)*by^2*dPsiBotSAy(i,:)./(2*kd(i)*PsiBotSAy(i,:));
    else
      GBotFunSAx(1,1,i,:) = sigma2*kplm^4*le*(lt+zbot-zbot.^3/(3*lb^2))./(4*(ks.^2-1/4*kd(i).^2));
      GBotFunSAy(1,1,i,:) = sigma2*kplm^4*le*(lt+zbot-zbot.^3/(3*lb^2))./(4*(ks.^2-1/4*kd(i).^2));
    end
end


FBotFunSAx = 1./sqrt(PsiBotSAx);
FBotFunSAy = 1./sqrt(PsiBotSAy);
FBotFunSA(1,1,:,:) = FBotFunSAx.*FBotFunSAy;


[ EXPtopSA(1,1,:,:)] = exp_topSA( sigma2,le,kd,kplm,ks,1/lt, ztop );
[ EXPbotSA(1,1,:,:)] = exp_botSA( sigma2,le,lb,kd,kplm,ks,1/lt, zbot );

PsiUnderSAx = ones(1,1,length(kd),length(zunder));
PsiUnderSAy = ones(1,1,length(kd),length(zunder));
dPsiUnderSAx = ones(1,1,length(kd),length(zunder));
dPsiUnderSAy = ones(1,1,length(kd),length(zunder));
GUnderFunSAx = ones(1,1,length(kd),length(zunder));
GUnderFunSAy = ones(1,1,length(kd),length(zunder));
EXPUnderSA = ones(1,1,length(kd),length(zunder));

for i = 1:length(kd)
    PsiUnderSAx(1,1,i,:) = PsiBotSAx(i,end)*(1+(zunder-lb)*dPsiBotSAx(i,end)./PsiBotSAx(i,end));
    PsiUnderSAy(1,1,i,:) = PsiBotSAy(i,end)*(1+(zunder-lb)*dPsiBotSAy(i,end)./PsiBotSAy(i,end));
    dPsiUnderSAx(1,1,i,:) = dPsiBotSAx(i,end);
    dPsiUnderSAy(1,1,i,:) = dPsiBotSAy(i,end);
    GUnderFunSAx(1,1,i,:) = GBotFunSAx(1,1,i,end)./(PsiUnderSAx(1,1,i,:)./PsiBotSAx(i,end));
    GUnderFunSAy(1,1,i,:) = GBotFunSAy(1,1,i,end)./(PsiUnderSAy(1,1,i,:)./PsiBotSAy(i,end));
    EXPUnderSA(1,1,i,:) =  exp(-sigma2*le*kd(i)^2*kplm^4/(8*(ks.^2-1/4*kd(i).^2).^2)*(2/3*lb+lt));
end

FUnderFunSAx = 1./sqrt(PsiUnderSAx);
FUnderFunSAy = 1./sqrt(PsiUnderSAy);
FUnderFunSA = FUnderFunSAx.*FUnderFunSAy;

GammaTop = ones(length(x),length(y),length(kd),length(ztop));
GammaBot = ones(length(x),length(y),length(kd),length(zbot));
GammaUnder = ones(length(x),length(y),length(kd),length(zunder));


for ix = 1:length(x)
    for iy = 1:length(y)
        GammaTop(ix,iy,:,:) = FTopFunSA(1,1,:,:).*exp(-x(ix).^2/ax^2*GTopFunSAx(1,1,:,:)-y(iy)^2/by^2*GTopFunSAy(1,1,:,:)).*EXPtopSA(1,1,:,:);
        GammaBot(ix,iy,:,:) = FBotFunSA(1,1,:,:).*exp(-x(ix).^2/ax^2*GBotFunSAx(1,1,:,:)-y(iy)^2/by^2*GBotFunSAy(1,1,:,:)).*EXPbotSA(1,1,:,:);
        GammaUnder(ix,iy,:,:) = FUnderFunSA(1,1,:,:).*exp(-x(ix).^2/ax^2*GUnderFunSAx(1,1,:,:)-y(iy)^2/by^2*GUnderFunSAy(1,1,:,:)).*EXPUnderSA(1,1,:,:);
    end
end
Coherency_Function = cat(4,GammaTop,GammaBot,GammaUnder);

name = ['Cartesian_Cyl_McDnld_phi' num2str(round(phi/pi*180)) '_AX' num2str(ax) '_BY' num2str(ax) '_sigma' num2str(sigma2) '_fs' num2str(fs/1000000000) 'HHz.mat'];
save(name)