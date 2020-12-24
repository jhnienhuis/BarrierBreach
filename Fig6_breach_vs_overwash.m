
addpath('D:\Dropbox\_Tools\delft3d_matlab')
load('PeakRoughness');

inputdir=[runsdir filesep 'ModelRun1' filesep];

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

idx = 27; %breach
idx = 4; %overwash

%adjust duration
    storm_surge = [1 output.storm_peak(idx) output.storm_peak(idx) 1 1];
    storm_time = [0 (12-output.duration(idx))*60 1+(12+output.duration(idx))*60 24*60 26*60];
    
    rundir=[runsdir filesep runname '_' num2str(idx,'%1.0f')];
    rundir='PeakRoughness_4';
    trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
    trih = vs_use([rundir filesep 'trih-bypass.dat'],[rundir filesep 'trih-bypass.def'],'quiet');
    t = vs_let(trim,'map-infsed-serie',{0},'MORFT','quiet');
    
    
    %roughness_map = vs_let(trim,'map-series',{0},'CFUROU',{0,0},'quiet');
    
    
    xc=vs_get(trim,'map-const','XCOR','quiet!');
    yc=vs_get(trim,'map-const','YCOR','quiet!');
    
    cell_area = cellarea(xc,yc);
    xcc = xc(:,1);
    ycc = yc(1,:);
    
    tt = round([0.01 0.25 0.5 0.75 1].*145);
    y_throat = 57;
    
    zmap = -1*vs_let(trim,'map-sed-series',{0},'DPS',{0,0},'quiet');
    
    z = -1*vs_let(trim,'map-sed-series',{tt},'DPS',{0,y_throat},'quiet')';
        
    h = vs_let(trim,'map-series',{tt},'S1',{0,y_throat},'quiet')';
    
    %calculate total volume of sediment
    qw = vs_let(trih,'his-series',{0},'CTR',{1},'quiet'); %m3s-1
    qs_bed = vs_let(trih,'his-sed-series',{0},'SBTR',{1,1},'quiet');
    qs_sus = vs_let(trih,'his-sed-series',{0},'SSTR',{1,1},'quiet'); %kgm-3
            
    dz = permute(zmap-zmap(1,:,:),[2 3 1]);
    dz(abs(dz)>10) = 0;
    idv_x = find(zmap(1,:,2)>1,1,'last')+(0:30); %170;
    idv_y = 57+[-4:4]; %2:111; %
    
    v = squeeze(sumdims(cell_area(idv_x,idv_y).*dz(idv_x,idv_y,:),[1 2]));
    zmap = zmap(tt,:,:);
    %plot(xcc,z), hold on, plot(xcc,h);

a = tight_subplot(2,5,0.01,0.01);
for ii=1:length(a), hold(a(ii),'on'), box(a(ii),'on'), end
cmap = demcmap([-2 3]);

h(:,end) = max(z(:,end),0);

for kk=1:length(tt)
    
    plot(a(kk),xcc,[z(:,1)],':k');
    plot(a(kk),xcc,[z(:,kk), h(:,kk)]);
    
    
    pcolor(a(kk+length(tt)),xc,yc,squeeze(zmap(kk,:,:))), shading(a(kk+length(tt)),'flat'), set(a(kk+length(tt)),'CLim',[-2 3])
    
end
for ii=1:length(a), set(a(ii),'XLim',[0 500]), end
for ii=(length(tt)+1):length(a), set(a(ii),'YLim',[-200 200],'DataAspectRatio',[1 1 1]), end
for ii=1:length(tt), set(a(ii),'YLim',[-3 3]), end
for ii=1:length(a)-1, set(a(ii),'PlotBoxAspectRatio',get(a(end),'PlotBoxAspectRatio')), end

for ii=1:length(a), set(a(ii),'Layer','top','XTickLabels','','YTickLabels',''),  end
colormap(cmap)
colorbar(a(7),'location','east')
saveas(gcf,'Fig6_breach_snaps.svg')

figure,
yyaxis('left')
plot(t*24,[qs_sus,qs_bed]), ylim([-0.01 0.05]),hold on
yyaxis('right')
plot(t*24,[qw v]), xlim([0 24]), ylim([-400 2000]), legend('water discharge','sus. sed transport','bed sed transport','overwash volume')
saveas(gcf,'Fig6_breach_timeseries.svg')