%% Сравнение функций когерентности для разных зависимостей электронной плотности и угла
%addpath('F:\Данные когерентности для статьи\Экватор\Результаты для ks и kd')
%load('phi101_AX25000_sigma001_fs1HHz.mat')
color = '-r';

disp('введите для какого значения x - переменной вывести данные')
disp('Допустимые значения ')
disp(y')
y1 = input('y = ');
ky=1;

disp('введите для какого значения fd - переменной вывести данные')
disp('Допустимые значения ')
disp(fd/1000000)
fd1 = input('fd = ');
fchoose=1;
   
   
while y(ky)~=y1
    ky=ky+1;   
end


while fd(fchoose)/1000000~=fd1
    fchoose=fchoose+1;   
end
  
module_x_01 = squeeze((abs(Coherency_Function(:,ky,fchoose,end))));
angle_x_01 = squeeze((angle(Coherency_Function(:,ky,fchoose,end))));
disp(fchoose)

fig1x = figure(1);
hold on
plot(x/1000,module_x_01,color,'LineWidth',2)
%title('Модуль функции когерентности SA') 
xlabel('x [km]','FontSize',16);
ylabel('|\Gamma|','FontSize',16)
set(gca,'Fontsize',16)
grid minor

fig2x = figure(2);
hold on
plot(x/1000,angle_x_01,color,'LineWidth',2)
%title('Модуль функции когерентности SA') 
xlabel('x [km]','FontSize',16);
ylabel('arg \Gamma','FontSize',16)
set(gca,'Fontsize',16)
grid minor

%   
  % q1 = legend(['$\varphi = ' num2str(74) '^\circ $ '],['$\varphi = ' num2str(86) '^\circ $ '],['$\varphi = ' num2str(101) '^\circ $ ']);
% q1 = legend('$\l_{\parallel}$ = 50 km','$\l_{\parallel}$ = 25 km','$\l_{\parallel}$ = 10 km');
% q1 = legend('$\sigma^2_N $ = ' num2str(sigma2) ,'$\sigma^2_N $ = ' num2str(sigma2) ,'$\sigma^2_N $ = ' num2str(sigma2) );
% set(q1,'interpreter','latex')
  
%     figure(2);
%     hold on
%     plot(x/1000,unwrap(angle11),'-b','LineWidth',2);
%     
%     xlabel('x [km]','FontSize',16);
%     ylabel('arg \Gamma','FontSize',16);
%     phi1 = round(phi/pi*180);
%     grid minor
% %   title('Фаза функции когерентности SA')
    disp(['phi = ' num2str(phi/pi*180)])
    disp(['lt = ' num2str(lt)])
    disp(['lb = ' num2str(lb)])
    disp(['kplm = ' num2str(kplm)])

  % hold on
  