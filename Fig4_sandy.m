%% add sandy data
%add better dune geometry?
sandy = xlsread('Nienhuis_BarrierBreach_tables.xlsx','Sandy');
a = logical(sandy(:,7));

%https://cera.coastalrisk.live/s/e2e0
overwash_volume = sandy(:,11)*0.3; %from eli's data

overwash_volume(overwash_volume<1) = -1;

%slopes, overwash volume, barrier width
subplot(2,2,1)
histogram(sandy(:,5)-sandy(:,6),[0:0.5:6])
xlabel('minimum barrier height (m)')

subplot(2,2,2)
histogram(sandy(:,4))
xlabel('barrier width (m)')

subplot(2,2,3)
histogram(sandy(:,8),[0:0.5:6])
xlabel('maximum water level difference (m)')

subplot(2,2,4)
histogram(overwash_volume)
xlabel('overwash volume (m3)')

saveas(gcf,'D:\Dropbox\2021 BarrierBreach JGR\Fig4_sandy.svg')
