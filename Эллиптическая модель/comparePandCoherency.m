clear
close all
load('CohEarth.mat');

Fs = (0.6:0.1:1.5)*10^9;
T0 = 10^-9;
fc = 0.9*10^9;

p02 = zeros(length(Fs), length(fd));
for i = 1:length(Fs)
    for j = 1:length(fd)
       p02(i,j) = T0^2/(2*pi)*exp(-T0^2/2*(fc-Fs(i))^2)*exp(-T0^2*fd(j)^2/4);
    end
end

figure(1)
grid minor
plot(fd, squeeze(abs(CoherenceEarth(4, 101, 101, :))))
hold on
plot(fd,p02(4,:)./max(p02(4,:)), 'r')
disp(['max(P^2)',num2str(max(p02(4,:)))])
disp(['max(Gamma)',num2str(max(squeeze(abs(CoherenceEarth(4, 101, 101, :)))))])