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
   x1 = input('y = ');
   kx=1;
   while y(kx)~=x1
    kx=kx+1;   
   end
   
  [X,Y]=meshgrid(fd/1000000,z_all/1000);
  
  figure
  module = squeeze((abs(Coherency_Function(kx,kx,:,:))));
  angle1 = squeeze((angle(Coherency_Function(kx,kx,:,:))));
  surf(X,Y,module');
 % title('Модуль частотной функции когерентности SA') 
  ylabel('z [km]','FontSize',16);
  xlabel('f_d [MHz]','FontSize',16);
  zlabel('|Г|','FontSize',16)
  colormap ('jet')
  set(gca,'Fontsize',20)
  view(135,45)
  
  figure
   surf(X,Y,unwrap(angle1,[],2)');
   ylabel('z [km]','FontSize',20);
   xlabel('f_d [MHz]','FontSize',20);
   zlabel('angle \Gamma','FontSize',20);
   colormap ('jet')
   set(gca,'Fontsize',20)
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
   ylabel('z [km]','FontSize',16);
   xlabel('y [km]','FontSize',16);
   zlabel('|\Gamma|','FontSize',16);
   colormap ('jet')
   
   set(gca,'Fontsize',20)
   view(135,45)
   
   figure
   surf(X,Y,unwrap(angle1,[],2)');
   ylabel('z [km]','FontSize',16);
   xlabel('y [km]','FontSize',16);
   zlabel('angle \Gamma','FontSize',16);
   colormap ('jet')
   set(gca,'Fontsize',20)
   view(135,45)
   %title('Фаза пространственной функции когерентности SA') 
end