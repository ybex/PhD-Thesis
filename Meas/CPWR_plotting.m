% CPWR plotting

%% global turn on

hfig = openfig('090459_Global_TurnOn');
d = findall(hfig, '-property', 'xdata');
xydatas = arrayfun(@(h) get(h, {'xdata','ydata', 'type'}), d, 'Uniform', 0);
close(gcf)

% ESR up
x=xydatas{1,1}{1,1};
y=xydatas{1,1}{1,2}*1e9;

h=figure(1);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(parula)


plot(x, y, 'Marker','none','LineStyle','-','LineWidth',3,'Color',c(15, :));

%xlim([43.29, 43.415]);
%ylim([0,1]);
xlabel('TT / TB / CC (V)')
ylabel('I (nA)')
set(gca,'FontSize',18, 'LineWidth', 3);
set(gca,'color','none')

%% pinch-off

data = load('Q:\spin-QED\CPWR_Data\2016\CQED_CWPR6GHz_v10\20160509\092642_Sweep_TB_vs_I .mat');
close all
x = data.data.plots{2}.xMat;
y =  data.data.plots{2}.yMat;
z = data.data.plots{2}.zMat*1e9;

h=figure(2);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(viridis);

pcolor(x, y, z)
shading flat;
c=colorbar;
c.Label.String  ='I (nA)';

xlim([1.5, 4]);
ylim([1 4])
xlabel('CC (V)')
ylabel('TT / TB (V)')
set(gca,'FontSize',18, 'LineWidth', 3);
set(c,'FontSize',18, 'LineWidth', 2);
set(gca,'color','none')

%% Resonance
% 0V
hfig = openfig('084257_Sweep-TB-I-0V');
d = findall(hfig, '-property', 'xdata');
xydatas = arrayfun(@(h) get(h, {'xdata','ydata', 'type'}), d, 'Uniform', 0);
close(gcf)

x=xydatas{1,1}{1,1}/1e9;
y=xydatas{1,1}{1,2};

% 4V
hfig = openfig('084335_Sweep-TB-I-4V');
d = findall(hfig, '-property', 'xdata');
xydatas = arrayfun(@(h) get(h, {'xdata','ydata', 'type'}), d, 'Uniform', 0);
close(gcf)

x2=xydatas{1,1}{1,1}/1e9;
y2=xydatas{1,1}{1,2};


h=figure(1);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(parula)

hold on
plot(x, y, 'Marker','none','LineStyle','-','LineWidth',3,'Color',c(5, :));
plot(x2, y2, 'Marker','none','LineStyle','-','LineWidth',3,'Color',c(35, :));

xlim([5.62, 5.7]);
ylim([-45,-30]);
xlabel('frequency (GHz)')
ylabel('Amplitude (dB)')
legend('CC = 0 V', 'CC = 4 V')
set(gca,'FontSize',18, 'LineWidth', 2);
set(gca,'color','none')

%% 2DEG sweep

data = load('090043_Sweep_TB_vs_I_Resonance .mat');
close all

% maximum
xm = data.data.plots{2}.xMat;
ym =  data.data.plots{2}.yMat;
zm = data.data.plots{2}.zMat;

% position
xp = data.data.plots{3}.xMat;
yp =  data.data.plots{3}.yMat;
zp = data.data.plots{3}.zMat;

% Q
xq = data.data.plots{4}.xMat;
yq =  data.data.plots{4}.yMat;
zq = data.data.plots{4}.zMat;

h=figure(2);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(viridis);

pcolor(xm, ym, zm)
shading flat;
c=colorbar;
c.Label.String  ='Resonance Amplitude (dB)';

xlim([1, 4]);
ylim([1 4])
xlabel('TT / TB (V)')
ylabel('CC (V)')
set(gca,'FontSize',18, 'LineWidth', 3);
set(c,'FontSize',18, 'LineWidth', 2);
set(gca,'color','none')

h=figure(3);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(viridis);

pcolor(xp, yp, zp)
shading flat;
c=colorbar;
c.Label.String  ='Resonance Position (GHz)';

xlim([1, 4]);
ylim([1 4])
xlabel('TT / TB (V)')
ylabel('CC (V)')
set(gca,'FontSize',18, 'LineWidth', 3);
set(c,'FontSize',18, 'LineWidth', 2);
set(gca,'color','none')

h=figure(4);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(viridis);

pcolor(xq, yq, zq)
shading flat;
c=colorbar;
%c.Label.String  ='Resonacne Q';

xlim([1, 4]);
ylim([1 4])
xlabel('TT / TB (V)')
ylabel('CC (V)')
set(gca,'FontSize',18, 'LineWidth', 3);
set(c,'FontSize',18, 'LineWidth', 2);
set(gca,'color','none')

%% Power 

data = load('191824_Sweep_P.mat');
close all
x = data.data.plots{2}.xMat/1e9;
y =  data.data.plots{2}.yMat;
z = data.data.plots{2}.zMat;

h=figure(2);
clf
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
c=colormap(viridis);

pcolor(x, y, z)
shading flat;
c=colorbar;
c.Label.String  ='Amplitude (dB)';

%xlim([1.5, 4]);
%ylim([1 4])
xlabel('frequency (GHz)')
ylabel('Power (dBm)')
set(gca,'FontSize',18, 'LineWidth', 3);
set(c,'FontSize',18, 'LineWidth', 2);
set(gca,'color','none')
