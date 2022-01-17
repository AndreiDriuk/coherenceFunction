% Построение функции когерентности
% addpath('F:\Данные когерентности для статьи\Экватор\Результаты для ks и kd')
% load('phi86_AX25000_sigma001_fs1HHz.mat')
close all
disp('Для построения частотной фукнции когерентности  введите 1' )
disp('Для построения пространственной функции когерентности  введите 2')
number=input('Number = ');
if number == 1
   disp('введите для какого значения x - переменной вывести данные')
   disp('Допустимые значения ')
   disp(x')
   x1 = input('x = ');
   kx=1;
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   disp('введите для какого значения y - переменной вывести данные')
   disp('Допустимые значения ')
   disp(y)
   y1 = input('y = ');
   ky=1;
   while y(ky)~=y1
    ky=ky+1;   
   end
   
  [X,Y]=meshgrid(fd/fs,z_all/1000);
  
  figure
  module = squeeze((abs(Coherency_Function(kx,ky,:,:))));
  angle1 = squeeze((angle(Coherency_Function(kx,ky,:,:))));
  surf(X,Y,module');
 % title('Модуль частотной функции когерентности SA') 
  ylabel('z [km]','FontSize',18);
  xlabel('f_d/f_s','FontSize',18);
  zlabel('|\Gamma|','FontSize',18)
  colormap ('jet')
  set(gca,'Fontsize',18)
  view(135,45)
  
  figure
   surf(X,Y,modunwrap(unwrap(angle1, [], 2)',3,2));
   %surf(X,Y,angle1');
   ylabel('z [km]','FontSize',18);
   xlabel(' f_d/f_s','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('Фаза частотной функции когерентности SA') 
   
elseif number==2
   disp('введите для какого значения частотной разностной переменной вывести данные')
   disp('Допустимые значения [МГЦ]')
   disp(fd/1000000)
   kd1=1;
   number1 = input('fd = ');
   while fd(kd1)/1000000~=number1
    kd1=kd1+1;   
   end
   
   disp('введите для какого значения y - переменной вывести данные')
   disp('Допустимые значения ')
   disp(y')
   y1 = input('y = ');
   ky=1;
   while y(ky)~=y1
    ky=ky+1;   
   end
   
   [X,Y]=meshgrid(x/1000,z_all/1000);
   
   figure
  module = squeeze((abs(Coherency_Function(:,ky,kd1,:))));
  angle1 = squeeze((angle(Coherency_Function(:,ky,kd1,:))));
  surf(X,Y,module');
   %title('Модуль пространственной функции когерентности SA') 
   ylabel('z [km]','FontSize',18);
   xlabel('x [km]','FontSize',18);
   zlabel('|\Gamma|','FontSize',18);
   colormap ('jet')
   
   set(gca,'Fontsize',20)
   view(135,45)
   
   figure
   surf(X,Y,modunwrap(unwrap(angle1,[],2)',3,2));
   ylabel('z [km]','FontSize',18);
   xlabel('x [km]','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('Фаза пространственной функции когерентности SA') 
   %%
   
   disp('введите для какого значения частотной разностной переменной вывести данные')
   disp('Допустимые значения [МГЦ]')
   disp(fd/1000000)
   kd1=1;
   number1 = input('fd = ');
   while fd(kd1)/1000000~=number1
    kd1=kd1+1;   
   end
   
   disp('введите для какого значения x - переменной вывести данные')
   disp('Допустимые значения ')
   disp(x')
   x1 = input('x = ');
   kx=1;
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   [X,Y]=meshgrid(y/1000,z_all/1000);
   
   figure
  module = squeeze((abs(Coherency_Function(kx,:,kd1,:))));
  angle1 = squeeze((angle(Coherency_Function(kx,:,kd1,:))));
  surf(X,Y,module');
   %title('Модуль пространственной функции когерентности SA') 
   ylabel('z [km]','FontSize',18);
   xlabel('y [km]','FontSize',18);
   zlabel('|\Gamma|','FontSize',18);
   colormap ('jet')
   
   set(gca,'Fontsize',20)
   view(135,45)
   
   figure
   surf(X,Y,modunwrap(unwrap(angle1,[],2)',3,2));
   ylabel('z [km]','FontSize',18);
   xlabel('y [km]','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('Фаза пространственной функции когерентности SA') 
   
end