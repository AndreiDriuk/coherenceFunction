Fs = 10^9;
Fd = ((-10:0.002:10)*10^8)';
lt = 140660; % coefficient in cosh(z/lt) 
lb = 198722; % coefficient in 1-z^2/lb^2
lend = 532000;
z = [-1000000:1000:lend,lend];
kplm =  0.2487;

charge = 1.602E-19;
vel = 299792458;
diel = 8.85E-12;
masse = 9.1E-31;
kplTop2 = kplm^2./cosh(z(z<=0)/lt);
kplBot2 = kplm^2*sqrt(1-z((z>0)&(z<=lb)).^2/lb^2);
kplUnder2 =z(z>lb)*0;
fpl2 = [kplTop2, kplBot2, kplUnder2]*vel^2/(4*pi^2);

phi1 = zeros(1, length(Fd));
phi2 = zeros(1, length(Fd));

tau_delay = trapz(z, (1./(1-fpl2/Fs^2).^(1/2) -1)/vel);

for j = 1:length(Fd)
        phi1(j) = (Fs+Fd(j)/2)/vel*trapz(z, (sqrt(1-fpl2/(Fs+Fd(j)/2).^2)-1));
        phi2(j) = (Fs-Fd(j)/2)/vel*trapz(z, (sqrt(1-fpl2/(Fs-Fd(j)/2).^2)-1));
end


%phimod1 = 2*pi/vel*Fd.*(1-fpl2/Fs^2)^(-3/4)*(1-fpl2/(Fs^2));

figure

plot(Fd,(abs(phi1-phi2-tau_delay.*Fd'))./abs(phi1-phi2))
set(gca, 'FontSize', 14)
grid minor
xlabel('f_d Гц')
ylabel('abs(phi1-phi2-tau_d*fd)/abs(phi1-phi2)')