e= 1.602176565e-19; %C 
h= 6.62606957e-34; %J*s
k_B = 1.38064852e-23; %m2kgs-2K-1
epsilon_0 = 8.854187817620e-12;
epsilon_r = 11.7;
Z = 1;

%% Donor potential for readout with electric field
h=figure(1)
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
c=colormap(parula)

N=10000;
r_start =0.1e-9;
r_stop = 25e-9;
r_pos = r_start:(r_stop-r_start)/N:r_stop;
r_tot = cat(2, -r_pos, r_pos);
pot = zeros(2*N+2,1); 
for i=1:2*N+2
    pot(i) = -Z*e^2/(4*pi*epsilon_0*epsilon_r*abs(r_tot(i)+1.5e-8))-5e6*r_tot(i)*e*1;
    pot_noE(i) = -Z*e^2/(4*pi*epsilon_0*epsilon_r*abs(r_tot(i)+1.5e-8));
    donor_E(i)=-5e6*(-1.5e-8)*1*e-0.045*e;
    donor_noE(i)=-0.045*e;
end
hold on
plot(r_tot*1e9, pot_noE/e*1000, 'Marker','.','LineStyle','none','LineWidth',2,'Color',c(30, :));
plot(r_tot*1e9, pot/e*1000, 'Marker','.','LineStyle','none','LineWidth',2,'Color',c(6, :));
plot(r_tot(4000:8000)*1e9, donor_E(4000:8000)/e*1000, 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(20, :));
plot(r_tot(4000:8000)*1e9, donor_noE(4000:8000)/e*1000, 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(45, :));
line([-25,20],[-30,-30], 'Marker','none','LineStyle','-.','LineWidth',2,'Color','k')
line([-25,20],[10,10], 'Marker','none','LineStyle','-.','LineWidth',2,'Color',c(1, :))
xlim([-25,0])
ylim([-220,120])
xlabel('r (nm)')
ylabel('E (meV)')
legend('E=0 MV/m', 'E=5 MV/m', 'D^0 E=0 MV/m', 'D^0 E=5 MV/m', 'Location', 'southeast')
set(gca,'FontSize',18, 'LineWidth', 3);
set(gca,'color','none')

%Schottky effect at the metal pulls down the conduction band if
%galvanically connected (yes, through GND)

%% gaussian for donor implantation
h=figure(2)
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
c=colormap(parula)
mu = 20;
sigma = 9.1;
skew = 0.49;
amplitude = 4e5;
x = 0:1:50;
skew = 2*normpdf(x,mu,sigma).*normcdf(skew*x)*amplitude;
figure(2)
plot(x, skew/1e4, 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(8, :))
set(gca,'FontSize',18, 'LineWidth', 3);
set(gca,'color','none')
xlabel('r (nm)')
ylabel('Atoms/\mu m^2')
