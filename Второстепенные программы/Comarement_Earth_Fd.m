%% Сравнение функций когерентности на Земле, как функций разностной частоты
%addpath('F:\Данные когерентности для статьи\Экватор\Результаты для ks и kd')
%load('phi86_AX25000_sigma01_fs1HHz.mat')
%color = '-b';

disp('введите для какого значения x - переменной вывести данные')
disp('Допустимые значения ')
disp(y')
y1 = input('y = ');
ky=1;

disp('введите для какого значения x - переменной вывести данные')
disp('Допустимые значения ')
disp(x')
x1 = input('x = ');
kx=1;


while y(ky)~=y1
    ky=ky+1;   
end


while x(kx)~=x1
    kx=kx+1;   
end

module_fd_01 = squeeze((abs(Coherency_Function(kx,ky,:,end))));
angle1_fd_01 = squeeze((angle(Coherency_Function(kx,ky,:,end))));
angle1_fd_01 = modunwrap(unwrap(angle1_fd_01), 2, 2);
% 
figure(5);
hold on
plot(fd/fs,module_fd_01, color,'LineWidth',2)
xlabel('f_d/f_s','FontSize',16);
ylabel('|Г|','FontSize',16)
set(gca,'Fontsize',16)
grid minor

figure(6);
hold on
plot(fd/fs ,unwrap(angle1_fd_01), color,'LineWidth',2)
xlabel('f_d/f_s','FontSize',16);
ylabel('arg \Gamma','FontSize',16)
set(gca,'Fontsize',16)
grid minor

% q1 = legend(['$\varphi = ' num2str(74) '^\circ $ '],['$\varphi = ' num2str(86) '^\circ $ '],['$\varphi = ' num2str(101) '^\circ $ ']);
% q1 = legend('$\l_{\parallel}$ = 50 km','$\l_{\parallel}$ = 25 km','$\l_{\parallel}$ = 10 km');
% q1 = legend('$\sigma^2_N $ = 0.0025','$\sigma^2_N $ = 0.01','$\sigma^2_N $ = 0.1');
% set(q1,'interpreter','latex')
 
% figure(2);
% hold on
% plot(fd/1000000,modunwrap(unwrap(angle11),3,2),'-b','LineWidth',2);
% xlabel('f_d [MHz]','FontSize',16);
% ylabel('arg \Gamma','FontSize',16);
% set(gca,'Fontsize',16)
% phi1 = round(phi/pi*180);

% q1 = legend(['$\varphi = ' num2str(74) '^\circ $ '],['$\varphi = ' num2str(86) '^\circ $ '],['$\varphi = ' num2str(101) '^\circ $ ']);
% q1 = legend('$\l_{\parallel}$ = 50 km','$\l_{\parallel}$ = 25 km','$\l_{\parallel}$ = 10 km');
% q1 = legend('$\sigma^2_N $ = 0.0025','$\sigma^2_N $ = 0.01','$\sigma^2_N $ = 0.1');
% set(q1,'interpreter','latex')

%set(gca,'Fontsize',16)
%grid minor
disp(['phi = ' num2str(phi/pi*180)])
disp(['lt = ' num2str(lt)])
disp(['lb = ' num2str(lb)])
disp(['kplm = ' num2str(kplm)])
disp(['AX = ' num2str(AX)] )
disp(['sigma2 = ' num2str(sigma2)] )

  