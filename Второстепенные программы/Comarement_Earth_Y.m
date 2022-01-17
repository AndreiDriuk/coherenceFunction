%% Сравнение функций когерентности для разных зависимостей электронной плотности и угла
%    addpath('F:\Данные когерентности для статьи\Экватор\Результаты для ks и kd')
%    load('phi86_AX25000_sigma01_fs1HHz.mat')

   % color = '-b';

   disp('введите для какого значения x - переменной вывести данные')
   disp('Допустимые значения ')
   disp(x')
   x1 = input('x = ');
   kx=1;
   
   disp('введите для какого значения fd - переменной вывести данные')
   disp('Допустимые значения ')
   disp(fd/1000000)
   fd1 = input('fd = ');
   fchoose=1;
   
   
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   
   while fd(fchoose)/1000000~=fd1
    fchoose=fchoose+1;   
   end
  
  module_y_01 = squeeze((abs(Coherency_Function(kx,:,fchoose,end))));
  angle_y_01 = squeeze((angle(Coherency_Function(kx,:,fchoose,end))));
  disp(fchoose)
  
  figure(3);
  hold on
  plot(y/1000,module_y_01,color,'LineWidth',2)
 %title('Модуль функции когерентности SA') 
  xlabel('y [km]','FontSize',16);
  ylabel('|Г|','FontSize',16)
  set(gca,'Fontsize',16)
  grid minor
  
  figure(4);
  hold on
  plot(y/1000,angle_y_01,color,'LineWidth',2)
 %title('Модуль функции когерентности SA') 
  xlabel('y [km]','FontSize',16);
  ylabel('arg \Gamma','FontSize',16)
  set(gca,'Fontsize',16)
  grid minor
% q1 = legend(['$\varphi = ' num2str(74) '^\circ $ '],['$\varphi = ' num2str(86) '^\circ $ '],['$\varphi = ' num2str(101) '^\circ $ ']);
% q1 = legend('$\l_{\parallel}$ = 50 km','$\l_{\parallel}$ = 25 km','$\l_{\parallel}$ = 10 km');
% q1 = legend('$\sigma^2_N $ = ' num2str(0.0025) ,'$\sigma^2_N $ = ' num2str(0.01) ,'$\sigma^2_N $ = ' num2str(0.1) );
% set(q1,'interpreter','latex')
  
%    figure(2);
%    hold on
%    plot(y/1000,unwrap(angle11),'-b','LineWidth',2);
%    
%    xlabel('y [km]','FontSize',16);
%    ylabel('arg \Gamma','FontSize',16);
%    phi1 = round(phi/pi*180);
%    grid minor
  %title('Фаза функции когерентности SA')
   disp(['phi = ' num2str(phi/pi*180)])
   disp(['lt = ' num2str(lt)])
   disp(['lb = ' num2str(lb)])
   disp(['kplm = ' num2str(kplm)])
   disp(['AX = ' num2str(AX)] )
   disp(['sigma2 = ' num2str(sigma2)] )
  