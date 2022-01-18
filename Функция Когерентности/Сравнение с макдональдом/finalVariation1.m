clear 
close all
addpath('../Дополнительные функции для расчета/')

ztop = -1000000:10000:0;
lt = 140660;
lb = 198722;
lend =  532000;
kplm  =  0.2487;
phi = 87.6274*pi/180;
zbot = [0:10000:lb,lb];
zunder = [lb:10000:lend,lend];
z_all = [ztop,zbot,zunder];
p = 11/3;
z1 = 0:50000;

fs=10^9;
fd = (0)';
sigma2=0.1;

x = -10000:1:10000;
y = -10000:1:10000;

v = 299792458; 
ks = fs*2*pi/v;
kd = fd*2*pi/v;
axm = 25000;
bym = 3000;
%lem = 671.6182;
AX = 9752;
BY = 144;
axp = sqrt((AX^2*BY^2)/((BY*sin(phi))^2+(AX*cos(phi))^2));
byp = BY;
lep = sqrt(AX^2*BY^2/((AX*sin(phi))^2+(BY*cos(phi))^2));

fun1 = lt*(1+tanh(ztop/lt));
fun2 = lt+zbot-zbot.^3/(3*lb^2);
fun3 = ones(1, length(zunder))*(lt+2/3*lb);
fun = [fun1, fun2, fun3];


AexReal = (1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/(2^(-(p-2)/2))).*exp(-(x*sin(phi)).^2/axm^2).*zMcDonald((p-2)/2, 2*pi*sqrt((x*cos(phi)).^2)/bym);
AeyReal = (1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/(2^(-(p-2)/2))).*zMcDonald((p-2)/2, 2*pi*y/bym);
AezReal = (1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/(2^(-(p-2)/2))).*exp(-(z1*cos(phi)).^2/axm^2).*zMcDonald((p-2)/2, 2*pi*(z1*sin(phi)/bym));
lem = trapz(z1, AezReal);
coeffReal = kplm^4*sigma2*lem/((ks^2-1/4*kd^2).^2);
DeXReal = 2*coeffReal*(AexReal(fix(length(x)/2+1))-AexReal);
DeYReal = 2*coeffReal*(AeyReal(fix(length(y)/2+1))-AeyReal);
CoherenceXReal = exp(-(ks^2-1/4*kd^2)/8*DeXReal'.*fun);
CoherenceYReal = exp(-(ks^2-1/4*kd^2)/8*DeYReal'.*fun);


SxReal = trapz(x,CoherenceXReal(:, end));
SyReal = trapz(y,CoherenceYReal(:, end));
SxPar = 0;
SyPar = 0;

while (SyPar<=SyReal)
    axp = sqrt(AX^2*BY^2/((BY*sin(phi))^2+(AX*cos(phi))^2));
    byp = BY;
    lep = sqrt(AX^2*BY^2/((AX*sin(phi))^2+(BY*cos(phi))^2));
    coeffPar = kplm^4*sigma2*lep/((ks^2-1/4*kd^2).^2);
    AeyPar = 1-y.^2/byp^2;
    DeYPar = 2*coeffPar*(AeyPar(fix(length(y)/2+1))-AeyPar);
    
    AexPar = 1-x.^2/(axp)^2;
    DeXPar = 2*coeffPar*(AexPar(fix(length(x)/2+1))-AexPar);
    
    CoherenceXPar = exp(-(ks^2-1/4*kd^2)/8*DeXPar'.*fun);
    CoherenceYPar = exp(-(ks^2-1/4*kd^2)/8*DeYPar'.*fun);
    SxPar = trapz(x,CoherenceXPar(:,end) );
    SyPar = trapz(y, CoherenceYPar(:,end));
    
    disp([' Y = ' num2str(BY),  ' Area Y = ', num2str(SyPar) , ' Area real Y ', num2str(SyReal), ' Area X = ', num2str(SxPar), ' Real Area ' num2str(SxReal)] )
    BY = BY+1;
end
BY = BY-1;

while (SxPar<=SxReal)
    axp = sqrt(AX^2*BY^2/((BY*sin(phi))^2+(AX*cos(phi))^2));
    byp = BY;
    lep = sqrt(AX^2*BY^2/((AX*sin(phi))^2+(BY*cos(phi))^2));
    coeffPar = kplm^4*sigma2*lep/((ks^2-1/4*kd^2).^2);
    AexPar = 1-x.^2/(axp)^2;
    DeXPar = 2*coeffPar*(AexPar(fix(length(x)/2+1))-AexPar);
    CoherenceXPar = exp(-(ks^2-1/4*kd^2)/8*DeXPar'.*fun);
    CoherenceYPar = exp(-(ks^2-1/4*kd^2)/8*DeYPar'.*fun);
    SxPar = trapz(x,CoherenceXPar(:,end) );
    SyPar = trapz(y, CoherenceYPar(:,end));
    disp([' X = ' num2str(AX), ' Area X = ', num2str(SxPar), ' Real Area» ' num2str(SxReal), ' Area Y = ', num2str(SyPar) , ' Area real Y', num2str(SyReal)] )
    AX = AX+1;
end
AX = AX-1;




disp([' X = ' num2str(AX), ' Area along X = ', num2str(SxPar), ' Real Area ' num2str(SxReal)] )

disp(['l_e = ', num2str(lep)])
disp([' X = ' num2str(AX)] );
disp([' Y = ' num2str(BY)]);


CoherenceXReal150kmAbove = CoherenceXReal(:,86);
CoherenceXPar150kmAbove = CoherenceXPar(:,86);
CoherenceXRealPierceP = CoherenceXReal(:,101);
CoherenceXParPierceP = CoherenceXPar(:,101);

CoherenceYReal150kmAbove = CoherenceYReal(:,86);
CoherenceYPar150kmAbove = CoherenceYPar(:,86);
CoherenceYRealPierceP = CoherenceYReal(:,101);
CoherenceYParPierceP = CoherenceYPar(:,101);

CoherenceXRealEarth = CoherenceXReal(:,end);
CoherenceYRealEarth = CoherenceYReal(:,end);
CoherenceXParEarth = CoherenceXPar(:,end);
CoherenceYParEarth = CoherenceYPar(:,end);

%% =======================================  
figure
plot(x/1000,CoherenceXRealEarth,'LineWidth', 2 )
hold on
plot(x/1000,CoherenceXParEarth,'LineWidth', 2 )
set(gca,'FontSize', 16)
xlabel('x km')
ylabel('Coherence function')
grid minor
xlim([-10,10])
legend('Power law model', 'Parabolic model')
set(gca,'FontSize', 16)
  
figure
plot(y/1000,CoherenceYRealEarth(:,end),'LineWidth', 2)
hold on
plot(y/1000,CoherenceYParEarth(:,end),'LineWidth', 2)
xlabel('y km')
ylabel('Coherence function')
grid minor
xlim([-1,1])
legend('Power law model', 'Parabolic model')
set(gca,'FontSize', 16)

%% ==================================
figure
plot(x/1000,CoherenceXRealPierceP,'LineWidth', 2 )
hold on
plot(x/1000,CoherenceXParPierceP,'LineWidth', 2 )
set(gca,'FontSize', 16)
xlabel('x km')
ylabel('Coherence function')
grid minor
xlim([-10,10])
legend('Power law model', 'Parabolic model')
set(gca,'FontSize', 16)
  
figure
plot(y/1000,CoherenceYRealPierceP,'LineWidth', 2)
hold on
plot(y/1000,CoherenceYParPierceP,'LineWidth', 2)
xlabel('y km')
ylabel('Coherence function')
grid minor
xlim([-1,1])
legend('Power law model', 'Parabolic model')
set(gca,'FontSize', 16)
%% ===================================================

figure
plot(x/1000,CoherenceXReal150kmAbove,'LineWidth', 2 )
hold on
plot(x/1000,CoherenceXPar150kmAbove,'LineWidth', 2 )
set(gca,'FontSize', 16)
xlabel('x km')
ylabel('Coherence function')
grid minor
xlim([-15,15])
legend('Power law model', 'Parabolic model')
set(gca,'FontSize', 16)
  
figure
plot(y/1000,CoherenceYReal150kmAbove,'LineWidth', 2)
hold on
plot(y/1000,CoherenceYPar150kmAbove,'LineWidth', 2)
xlabel('y km')
ylabel('Coherence function')
grid minor
xlim([-1,1])
legend('Power law model', 'Parabolic model')
set(gca,'FontSize', 16)