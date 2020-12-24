
addpath('D:\Dropbox\_Tools\delft3d_matlab')
load('PeakRoughness');
inputdir=[runsdir filesep 'ModelRun1' filesep];

idx = 4;
rundir=[runsdir filesep runname '_' num2str(idx,'%1.0f')];
rundir='PeakRoughness_4';


%adjust duration
storm_surge = [1 output.storm_peak(idx) output.storm_peak(idx) 1 1];
storm_time = [0 (12-output.duration(idx))*60 1+(12+output.duration(idx))*60 24*60 26*60];


%note: these functions come with the delft3d installation at
%C:\Program Files (x86)\Deltares\Delft3D 4.01.00\win32\delft3d_matlab
trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
trih = vs_use([rundir filesep 'trih-bypass.dat'],[rundir filesep 'trih-bypass.def'],'quiet');
t = vs_let(trim,'map-infsed-serie',{0},'MORFT','quiet');

roughness_map = squeeze(vs_let(trim,'map-series',{50},'CFUROU',{0,0},'quiet'));
%imagesc(squeeze(roughness_map(80,:,:)),[10 70])


xc=vs_get(trim,'map-const','XCOR','quiet!');
yc=vs_get(trim,'map-const','YCOR','quiet!');

cell_area = cellarea(xc,yc);
xcc = xc(:,1);
ycc = yc(1,:);

zmap = squeeze(-1*vs_let(trim,'map-sed-series',{1},'DPS',{0,0},'quiet'))';


subplot(2,2,1)

pcolor(xcc,ycc(2:end)',zmap(2:end,:)), shading flat
colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
title(gca,'bathymetry')
xlim([0 2000])
colorbar('location','east')
box on
subplot(2,2,2)
pcolor(xcc(2:end),ycc(2:end-1)',roughness_map(2:end,2:end-1)'), shading flat, set(gca,'CLim',[10 50])
title(gca,'roughness')
xlim([0 2000])
box on
subplot(2,2,3)
pcolor(xcc,ycc,zeros(size(zmap)))
title(gca,'model grid')
xlim([0 2000])

subplot(2,2,4)
yyaxis('left')
plot(storm_time./60,storm_surge-1), xlim([0 24]), set(gca,'XTick',[0 6 12 18 24])
ylabel('ocean water level (m)')
yyaxis('right')
plot(storm_time./60,(storm_surge-1)./2000), xlim([0 24]), set(gca,'XTick',[0 6 12 18 24])
title(gca,'model grid')
xlabel('time (hr)')
ylabel('water level slope across domain')
saveas(gcf,'Fig3_modelsetup.svg')