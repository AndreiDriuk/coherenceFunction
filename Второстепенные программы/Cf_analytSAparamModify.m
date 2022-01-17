clear 
tic
ztop =[-100000,0];% -1000000:20000:0;
lt = 154154;
lb = 147536;
lend = 422000;
kplm = 0.20544;
phi = 86.0433/180*pi;
zbot = [0,lb];%[0:10000:lb,lb];
zunder = [lb,lend];%[lb:10000:lend,lend];
z_all = [ztop,zbot,zunder];

fs=10^9;
fd =([-6:0.5:6]*10^8)';
sigma2=0.01;

v = 299792458; %скорость распространения
ks = fs*2*pi/v;
kd = fd*2*pi/v;
% пространственный масштаб неоднородностей
x = -2000:100:2000;
y = -1000:50:1000;
leX = 25000;
leY = 3000;
ler = sqrt(leX^2*leY^2/((leX*sin(phi))^2+(leY*cos(phi))^2));%sqrt((leX*cos(phi))^2+(leY*sin(phi))^2);
c = leX/leY;
p = 11/3;

%% Определение эффективных масштабов
Ae = zeros(length(x), length(y));
Par = zeros(length(x), length(y));
for iy = 1:length(y)
    Ae(:,iy) =(1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/(2^(-(p-2)/2))).*exp(-(x*sin(phi)).^2/leX^2).*zMcDonald((p-2)/2, 2*pi*sqrt((x*cos(phi)).^2+y(iy)^2)/leY);
    Par(:,iy) = x.^2/(c^2/(c^2*cos(phi)^2+sin(phi)^2))+y(iy)^2;
end
De = (Ae(fix(length(x)/2+1),fix(length(y)/2+1))-Ae);

BY = sqrt(c^2/((c*sin(phi))^2+cos(phi)^2))*Par./(ler.*De);
BY((x==0),(y==0)) = ler;
AX = c*BY;
ax = c*BY*sqrt((1/(c^2*cos(phi)^2+sin(phi)^2)));
by = BY;
le = BY*sqrt(c^2/((c*sin(phi))^2+cos(phi)^2));
%%

nu_tSAx = zeros(length(x), length(y), length(kd));
nu_tSAy = zeros(length(x), length(y), length(kd));
nu_bx = zeros(length(x), length(y), length(kd));
nu_by = zeros(length(x), length(y), length(kd));
ksi_bx = zeros(length(x), length(y), length(kd), length(zbot));
ksi_by = zeros(length(x), length(y), length(kd), length(zbot));
dksi_bx = zeros(length(x), length(y), length(kd), length(zbot));
dksi_by = zeros(length(x), length(y), length(kd), length(zbot));

for ix=1:length(x)
    for iy = 1:length(y)
        nu_tSAx(ix, iy, :)=-1/2-1/2*sqrt(1+2*1i*kd*sigma2*kplm^4*le(ix, iy)./((ks.^2-1/4*kd.^2).^2/lt^2*ax(ix, iy)^2));
        nu_tSAy(ix, iy, :)=-1/2-1/2*sqrt(1+2*1i*kd*sigma2*kplm^4*le(ix, iy)./((ks.^2-1/4*kd.^2).^2/lt^2*by(ix, iy)^2));
        nu_bx(ix, iy, :)=sqrt(sigma2)*kplm^2*lb*sqrt(1i*kd*le(ix, iy)/(2*ax(ix, iy)^2))./(2*(ks.^2-1/4*kd.^2))-1/2;
        nu_by(ix, iy, :)=sqrt(sigma2)*kplm^2*lb*sqrt(1i*kd*le(ix, iy)/(2*by(ix, iy)^2))./(2*(ks.^2-1/4*kd.^2))-1/2;
        for ifd = 1:length(fd)
            ksi_bx(ix,iy,ifd,:)=zbot.*sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(ifd)*le(ix, iy)/(2*ax(ix, iy)^2))./((ks.^2-1/4*kd(ifd).^2)*lb));
            ksi_by(ix,iy,ifd,:)=zbot.*sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(ifd)*le(ix, iy)/(2*by(ix, iy)^2))./((ks.^2-1/4*kd(ifd).^2)*lb));
            dksi_bx(ix,iy,ifd,:)=sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(ifd)*le(ix, iy)/(2*ax(ix, iy)^2))./((ks.^2-1/4*kd(ifd).^2)*lb));
            dksi_by(ix,iy,ifd,:)=sqrt(2*sqrt(sigma2)*kplm^2*sqrt(1i*kd(ifd)*le(ix, iy)/(2*by(ix, iy)^2))./((ks.^2-1/4*kd(ifd).^2)*lb));
        end
    end
end


PsiTopSAx = zeros(length(x), length(y), length(kd), length(ztop));
dPsiTopSAx = zeros(length(x), length(y), length(kd), length(ztop));
PsiTopSAy = zeros(length(x), length(y), length(kd), length(ztop));
dPsiTopSAy = zeros(length(x), length(y), length(kd), length(ztop));
for ix = 1:length(x)
    for iy = 1:length(y)
        [PsiTopSAx(ix, iy,:,:), dPsiTopSAx(ix, iy,:,:)] = fhitop(squeeze(nu_tSAx(ix, iy,:,:)),1/lt,ztop);
        [PsiTopSAy(ix, iy,:,:), dPsiTopSAy(ix, iy,:,:)] = fhitop(squeeze(nu_tSAy(ix, iy,:,:)),1/lt,ztop);
    end
end


GTopFunSAx = zeros(length(x), length(y),length(kd),length(ztop));
GTopFunSAy = zeros(length(x), length(y),length(kd),length(ztop));
for ix = 1:length(x)
    for iy = 1:length(y)
        for i =1:length(kd)
            if kd(i)~=0
                GTopFunSAx(ix,iy,i,:) = 1i*ax(ix, iy)^2*(ks.^2-1/4*kd(i).^2)*dPsiTopSAx(ix,iy,i,:)./(2*kd(i)*PsiTopSAx(ix,iy,i,:));
                GTopFunSAy(ix,iy,i,:) = 1i*by(ix, iy)^2*(ks.^2-1/4*kd(i).^2)*dPsiTopSAy(ix,iy,i,:)./(2*kd(i)*PsiTopSAy(ix,iy,i,:));
            else
                hyper=hypergeom([squeeze(-nu_tSAx(ix, iy, i)+1),squeeze(nu_tSAx(ix, iy,i)+2)],2,((1+tanh(ztop/lt))/2));
                GTopFunSAx(ix,iy,i,:) = sigma2*kplm^4*le(ix, iy).*hyper./(8*(ks.^2-1/4*kd(i).^2)/lt*cosh(ztop/lt).^2.*squeeze(PsiTopSAx(ix,iy,i,:)).');
                hyper=hypergeom([-nu_tSAy(ix, iy,i)+1,nu_tSAy(ix, iy,i)+2],2,((1+tanh(ztop/lt))/2));
                GTopFunSAy(ix,iy,i,:) = sigma2*kplm^4*le(ix, iy).*hyper./(8*(ks.^2-1/4*kd(i).^2)/lt*cosh(ztop/lt).^2.*squeeze(PsiTopSAy(ix,iy,i,:)).');      
            end
        end
    end
end
 clear hyper nu_tSAx nu_tSAy


FTopFunSAx=1./sqrt(PsiTopSAx);
FTopFunSAy=1./sqrt(PsiTopSAy);
FTopFunSA=FTopFunSAx.*FTopFunSAy;

DEV_SAx = zeros(length(x), length(y),length(kd),length(zbot)); 
DEV_SAy = zeros(length(x), length(y),length(kd),length(zbot));
DOV_SAx = zeros(length(x), length(y),length(kd),length(zbot));
DOV_SAy = zeros(length(x), length(y),length(kd),length(zbot));
DDEV_SAx = zeros(length(x), length(y),length(kd),length(zbot));
DDEV_SAy = zeros(length(x), length(y),length(kd),length(zbot));
DDOV_SAx = zeros(length(x), length(y),length(kd),length(zbot));
DDOV_SAy = zeros(length(x), length(y),length(kd),length(zbot));
PsiBotSAx = zeros(length(x), length(y),length(kd),length(zbot));
PsiBotSAy = zeros(length(x), length(y),length(kd),length(zbot));
dPsiBotSAx = zeros(length(x), length(y),length(kd),length(zbot));
dPsiBotSAy = zeros(length(x), length(y),length(kd),length(zbot));
for ix = 1:length(x)
    for iy = 1:length(y)
        [ DEV_SAx(ix,iy,:,:), DOV_SAx(ix,iy,:,:),DDEV_SAx(ix,iy,:,:),DDOV_SAx(ix,iy,:,:)] = weberbot(squeeze(nu_bx(ix,iy,:)),squeeze(ksi_bx(ix,iy,:,:)),squeeze(dksi_bx(ix,iy,:,:)));
        [ DEV_SAy(ix,iy,:,:), DOV_SAy(ix,iy,:,:),DDEV_SAy(ix,iy,:,:),DDOV_SAy(ix,iy,:,:)] = weberbot(squeeze(nu_by(ix,iy,:)),squeeze(ksi_by(ix,iy,:,:)),squeeze(dksi_by(ix,iy,:,:)));
        [ PsiBotSAx(ix,iy,:,:), dPsiBotSAx(ix,iy,:,:) ] = phibot(squeeze(DEV_SAx(ix,iy,:,:)),squeeze(DOV_SAx(ix,iy,:,:)),squeeze(DDEV_SAx(ix,iy,:,:)),squeeze(DDOV_SAx(ix,iy,:,:)),squeeze(PsiTopSAx(ix,iy,:,:)),squeeze(dPsiTopSAx(ix,iy,:,:)),kd,zbot);
        [ PsiBotSAy(ix,iy,:,:), dPsiBotSAy(ix,iy,:,:) ] = phibot(squeeze(DEV_SAy(ix,iy,:,:)),squeeze(DOV_SAy(ix,iy,:,:)),squeeze(DDEV_SAy(ix,iy,:,:)),squeeze(DDOV_SAy(ix,iy,:,:)),squeeze(PsiTopSAy(ix,iy,:,:)),squeeze(dPsiTopSAy(ix,iy,:,:)),kd,zbot);
    end
end
% 
% 
 clear DEV_SAx DOV_SAx DDEV_SAx DDOV_SAx DEV_SAy DOV_SAy DDEV_SAy DDOV_SAy
 GBotFunSAx = ones(length(x),length(y),length(kd),length(zbot));
 GBotFunSAy = ones(length(x),length(y),length(kd),length(zbot));
for ix = 1:length(x)
    for iy = 1:length(y)
        for i =1:length(kd)
            if kd(i)~=0
                GBotFunSAx(ix,iy,i,:) = 1i*(ks.^2-1/4*kd(i).^2)*ax(ix, iy)^2*dPsiBotSAx(ix, iy,i,:)./(2*kd(i)*PsiBotSAx(ix, iy,i,:));
                GBotFunSAy(ix,iy,i,:) = 1i*(ks.^2-1/4*kd(i).^2)*by(ix, iy)^2*dPsiBotSAy(ix, iy,i,:)./(2*kd(i)*PsiBotSAy(ix, iy,i,:));
            else
                GBotFunSAx(ix,iy,i,:) = sigma2*kplm^4*le(ix, iy)*(lt+zbot-zbot.^3/(3*lb^2))./(4*(ks.^2-1/4*kd(i).^2));
                GBotFunSAy(ix,iy,i,:) = sigma2*kplm^4*le(ix, iy)*(lt+zbot-zbot.^3/(3*lb^2))./(4*(ks.^2-1/4*kd(i).^2));
            end
        end
    end
end

 
FBotFunSAx = 1./sqrt(PsiBotSAx);
FBotFunSAy = 1./sqrt(PsiBotSAy);
FBotFunSA = FBotFunSAx.*FBotFunSAy;
EXPtopSA = ones(length(x),length(y),length(kd),length(ztop));
EXPbotSA = ones(length(x),length(y),length(kd),length(zbot));
for ix = 1:length(x)
    for iy = 1:length(y)
        [ EXPtopSA(ix,iy,:,:)] = exp_topSA( sigma2,le(ix,iy),kd,kplm,ks,1/lt, ztop );
        [ EXPbotSA(ix,iy,:,:)] = exp_botSA( sigma2,le(ix,iy),lb,kd,kplm,ks,1/lt, zbot );
    end
end

 
PsiUnderSAx = ones(length(x),length(y),length(kd),length(zunder));
PsiUnderSAy = ones(length(x),length(y),length(kd),length(zunder));
dPsiUnderSAx = ones(length(x),length(y),length(kd),length(zunder));
dPsiUnderSAy = ones(length(x),length(y),length(kd),length(zunder));
GUnderFunSAx = ones(length(x),length(y),length(kd),length(zunder));
GUnderFunSAy = ones(length(x),length(y),length(kd),length(zunder));
EXPUnderSA = ones(length(x),length(y),length(kd),length(zunder));

for ix = 1:length(x)
    for iy = 1:length(y)
        for i = 1:length(kd)
            PsiUnderSAx(ix,iy,i,:) = PsiBotSAx(ix,iy,i,end)*(1+(zunder-lb)*dPsiBotSAx(ix,iy,i,end)./PsiBotSAx(ix,iy,i,end));
            PsiUnderSAy(ix,iy,i,:) = PsiBotSAy(ix,iy,i,end)*(1+(zunder-lb)*dPsiBotSAy(ix,iy,i,end)./PsiBotSAy(ix,iy,i,end));
            dPsiUnderSAx(ix,iy,i,:) = dPsiBotSAx(ix,iy,i,end);
            dPsiUnderSAy(ix,iy,i,:) = dPsiBotSAy(ix,iy,i,end);
            GUnderFunSAx(ix,iy,i,:) = GBotFunSAx(ix,iy,i,end)./(PsiUnderSAx(ix,iy,i,:)./PsiBotSAx(ix,iy,i,end));
            GUnderFunSAy(ix,iy,i,:) = (squeeze(GBotFunSAy(ix,iy,i,end)).')./((squeeze(PsiUnderSAy(ix,iy,i,:)).')./squeeze(PsiBotSAy(ix,iy,i,end)).');
            EXPUnderSA(ix,iy,i,:) =  exp(-sigma2*le(ix,iy)*kd(i)^2*kplm^4/(8*(ks.^2-1/4*kd(i).^2).^2)*(2/3*lb+lt));
        end
    end
end

FUnderFunSAx = 1./sqrt(PsiUnderSAx);
FUnderFunSAy = 1./sqrt(PsiUnderSAy);
FUnderFunSA = FUnderFunSAx.*FUnderFunSAy;

GammaTop = ones(length(x),length(y),length(kd),length(ztop));
GammaBot = ones(length(x),length(y),length(kd),length(zbot));
GammaUnder = ones(length(x),length(y),length(kd),length(zunder));


for ix = 1:length(x)
    for iy = 1:length(y)
        GammaTop(ix,iy,:,:) = FTopFunSA(ix,iy,:,:).*exp(-x(ix).^2/ax(ix,iy)^2*GTopFunSAx(ix,iy,:,:)-y(iy)^2/by(ix,iy)^2*GTopFunSAy(ix,iy,:,:)).*EXPtopSA(ix,iy,:,:);
        GammaBot(ix,iy,:,:) = FBotFunSA(ix,iy,:,:).*exp(-x(ix).^2/ax(ix,iy)^2*GBotFunSAx(ix,iy,:,:)-y(iy)^2/by(ix,iy)^2*GBotFunSAy(ix,iy,:,:)).*EXPbotSA(ix,iy,:,:);
        GammaUnder(ix,iy,:,:) = FUnderFunSA(ix,iy,:,:).*exp(-x(ix).^2/ax(ix,iy)^2*GUnderFunSAx(ix,iy,:,:)-y(iy)^2/by(ix,iy)^2*GUnderFunSAy(ix,iy,:,:)).*EXPUnderSA(ix,iy,:,:);
    end
end
Coherency_Function = cat(4,GammaTop,GammaBot,GammaUnder);
toc
% 
% name = ['McDnld_phi' num2str(round(phi/pi*180)) '_AX' num2str(AX) '_sigma' num2str(sigma2) '_fs' num2str(fs/1000000000) 'HHz.mat'];
% save(name)