%% add sandy data
%add better dune geometry?
sandy = xlsread('Nienhuis_BarrierBreach_tables.xlsx','Sandy');
a = logical(sandy(:,7));
d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

s = size(sandy,1);

for ii=1:s,
storm_time = sandy(ii,9)*3600.*[0 0.5 1 1.5 2];
storm_surge = [0 sandy(ii,8) sandy(ii,8) 0 0];

Cdrag = gravity./((55-(sandy(ii,7).*30)).^2);

E = predict_volume(storm_time,storm_surge,sandy(ii,8),sandy(ii,5)-3,Cdrag,gravity,Rsed,sandy(ii,4),d50);

Vs(ii) = E.*100*sandy(ii,4);
%Eobs = max(z(22,2:end-1,end))
end

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
saveas(gcf,'Fig4_sandy.svg')
