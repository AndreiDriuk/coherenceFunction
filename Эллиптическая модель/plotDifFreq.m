load('0--30BY50000sig001.mat')
   disp('введите для какого значения x - переменной вывести данные')
   disp('Допустимые значения ')
   disp(x')
   x1 = input('x = ');
   kx=1;

   
   disp('введите для какого значения y - переменной вывести данные')
   disp('Допустимые значения ')
   disp(y')
   y1 = input('y = ');
   ky=1;
   
   disp('введите для какого значения fd - переменной вывести данные')
   disp('Допустимые значения ')
   disp(fd/1000000)
   fd1 = input('fd = ');
   fchoose=1;
   
   
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   
   while y(ky)~=y1
    ky=ky+1;   
   end
   
   while fd(fchoose)/1000000~=fd1
    fchoose=fchoose+1;   
   end
  
  module1 = squeeze((abs(Coherency_Function(kx,ky,fchoose,:))));
  angle11 = squeeze((angle(Coherency_Function(kx,ky,fchoose,:))));
  disp(fchoose)
  
  fig1 = figure(1);
  hold on
  plot(z_all/1000,module1,'-b','LineWidth',2)
 %title('Модуль функции когерентности SA') 
  xlabel('z [km]','FontSize',16);
  ylabel('|Г|','FontSize',16)
  set(gca,'Fontsize',20)
  
   figure(2);
   hold on
   plot(z_all/1000,unwrap(angle11),'-b','LineWidth',2);
   
   xlabel('z [km]','FontSize',16);
   ylabel('phase \Gamma','FontSize',16);
   phi1 = round(phi/pi*180);
  %title('Фаза функции когерентности SA')
   disp(['phi = ' num2str(phi/pi*180)])
   disp(['lt = ' num2str(lt)])
   disp(['lb = ' num2str(lb)])
   disp(['kplm = ' num2str(kplm)])

   load('0--30BY50000sig001f108.mat')
   disp(['phi = ' num2str(phi/pi*180)])
   disp(['lt = ' num2str(lt)])
   disp(['lb = ' num2str(lb)])
   disp(['kplm = ' num2str(kplm)])
  module1 = squeeze((abs(Coherency_Function(kx,ky,fchoose,:))));
  angle11 = squeeze((angle(Coherency_Function(kx,ky,fchoose,:))));
     phi2 = round(phi/pi*180);
  disp(fchoose)
  
    figure(1);
    hold on
   plot(z_all/1000,module1,'--r','LineWidth',2)
  
   figure(2);
   hold on
   plot(z_all/1000,(unwrap(angle11))','--r','LineWidth',2);
   
   
   figure(1);
   hold on
   plot(z_all/1000,module1,'-.g','LineWidth',2)
   grid on
   set(gca,'FontSize',14)
   q = legend('f = 1 HHz','f = 100 MHz');
   
   figure(2);
   grid on
   set(gca,'FontSize',14)
   q1 = legend('f = 1 HHz','f = 100 MHz');

   