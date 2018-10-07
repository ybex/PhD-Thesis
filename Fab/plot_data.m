HF4K = csvread('HF4K.csv');
HFRT = csvread('HFRT.csv');
LF4K = csvread('LF4K.csv');
LFRT = csvread('LFRT.csv');

%%
h=figure(203);
clf
hold all
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
c=colormap(parula)

plot(HF4K(:,1)/1e9,HF4K(:,2), 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(1, :))
plot(HFRT(:,1)/1e9,HFRT(:,2), 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(20, :))
plot(LF4K(:,1)/1e9,LF4K(:,2), 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(38,:))
plot(LFRT(:,1)/1e9,LFRT(:,2), 'Marker','none','LineStyle','-','LineWidth',2,'Color',c(48,:))
legend('HF in, 4K','HF in, RT','LF in, 4K','LF in, RT', 'Location', 'south')
xlabel('Frequency (GHz)')
ylabel('S_{12} (dB)')
xlim([0 26.5])


set(gca,'FontSize',18, 'LineWidth', 3);

%legend('Nemo3D similation', 'theory')
set(gca,'color','none')

print('Diplexer18GHz.pdf','-r0','-dpdf')