%% breach movie
addpath('D:\Dropbox\_Tools\delft3d_matlab')
load('PeakRoughness');

inputdir=[runsdir filesep 'ModelRun1' filesep];

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

idx = 27; %breach
idx = 4; %overwash test

  

storm_surge = [1 output.storm_peak(idx) output.storm_peak(idx) 1 1];
storm_time = [0 (12-output.duration(idx))*60 1+(12+output.duration(idx))*60 24*60 26*60];

rundir=[runsdir filesep runname '_' num2str(idx,'%1.0f')];
savename = ['D:\Dropbox\2021 BarrierBreach JGR\' runname '_' num2str(idx,'%1.0f') '.gif'];

trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
trih = vs_use([rundir filesep 'trih-bypass.dat'],[rundir filesep 'trih-bypass.def'],'quiet');
t = vs_let(trim,'map-infsed-serie',{0},'MORFT','quiet');

xc=vs_get(trim,'map-const','XCOR','quiet!');
yc=vs_get(trim,'map-const','YCOR','quiet!');

cell_area = cellarea(xc,yc);
xcc = xc(:,1);
ycc = yc(1,:);

y_throat = 57;

z = permute(-1*vs_let(trim,'map-sed-series',{0},'DPS',{0,0},'quiet'),[2 3 1]);
h = vs_let(trim,'map-series',{0},'S1',{0,0},'quiet');

dz = z-z(:,:,1);
dz(abs(dz)>10) = 0;
idv_x = find(z(:,2,1)>1,1,'last')+(0:30); %170;
idv_y = 57+[-4:4]; %2:111; %

v = squeeze(sumdims(cell_area(idv_x,idv_y).*dz(idv_x,idv_y,:),[1 2]));

idv_x = find(z(:,y_throat,1)>0);
idv_y = find(z(20,:,1)<(output.height(idx)+0.1)); idv_y(idv_y==1 | idv_y==112) = [];

vbar = output.vbarrier(idx) + squeeze(sumdims(cell_area(idv_x,idv_y).*dz(idv_x,idv_y,:),[1 2]));


fig = figure('color','white','Name','Breach_WaterSurface','visible','off','Units','Pixels','PaperUnits','Centimeters','Units','Centimeter',...
    'InvertHardCopy', 'off','NextPlot','add','Position',[5 5 35 15]);

a = tight_subplot(1,2,0.02,0.02);
%a(1) = subplot(1,2,1);
%a(2) = subplot(1,2,2);

for ii=1:length(a), hold(a(ii),'on'), box(a(ii),'on'), end

cmap = demcmap([-3 3]);

%fixed plot
qs_cum = cumsum(1/(1-0.4).*output.qs_cum(idx,:))';
h2 = plot(a(2),t*24,[vbar qs_cum]);
plot(a(2),[0 24],[0 0],':k')
c = get(h2,'Color'); 


for kk=1:length(t)
    zmap = squeeze(z(:,:,kk));
    hmap = squeeze(h(kk,:,:));   
    nanmap = (-0.2+hmap)<=zmap;
    hmap = min(-0.1,1+(hmap.*-1));
    hmap(nanmap) = nan;
          
    h0 = pcolor(a(1),xc,yc,zmap);
    h1 = pcolor(a(1),xc,yc,hmap);
    
    shading(a(1),'flat'),
    set(a(1),'CLim',[-3 3]);
    set(a(1),'XLim',[0 500],'YLim',[-200 200],'DataAspectRatio',[1 1 1])
    set(a(1),'Layer','top','XTickLabels','','YTickLabels','')
    
    

    h3 = scatter(a(2),t(kk)*24,[vbar(kk) qs_cum(kk)],50,cell2mat(c),'filled');
    
    set(a(2),'Layer','top','xlim',[0 24],'YAxisLocation','right')
    xticks(a(2),[0 6 12 18 24])
    xticklabels(a(2),[0 6 12 18 24])
    yticklabels(a(2),'auto')
    set(a(2),'PlotBoxAspectRatio',get(a(1),'PlotBoxAspectRatio'))
    xlabel(a(2),'Time (hrs)')
    legend(a(2),'Barrier Volume (m^3)','Washover Flux (m^3)') 
    
    colormap(cmap)
    
    cb = colorbar(a(1),'location','east');
    cb.TickLabels = {'3','2','1','0','1','2','3'};
    xlabel(cb,'<- Water elevation (m) | Barrier elevation (m) ->')
    
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),64);
    if kk == 1;
        imwrite(imind,cm,savename,'gif', 'Loopcount',inf,'DelayTime',1/20);
    else
        imwrite(imind,cm,savename,'gif','WriteMode','append','DelayTime',1/20);
    end
    
    delete(h0)
    delete(h1)
    delete(h3)
    
    
end

imwrite(imind,cm,savename,'gif','WriteMode','append','DelayTime',1/5);



%% side by side breaching and overwash
addpath('D:\Dropbox\_Tools\delft3d_matlab')
load('PeakRoughness');

inputdir=[runsdir filesep 'ModelRun1' filesep];

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

idx = 4; %overwash
idx2 = 27; %breach

savename = ['D:\Dropbox\2021 BarrierBreach JGR\' runname '_' num2str(idx,'%1.0f') '_' num2str(idx2,'%1.0f') '.gif'];

rundir=[runsdir filesep runname '_' num2str(idx,'%1.0f')];
trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
t = vs_let(trim,'map-infsed-serie',{0},'MORFT','quiet');

xc=vs_get(trim,'map-const','XCOR','quiet!');
yc=vs_get(trim,'map-const','YCOR','quiet!');

cell_area = cellarea(xc,yc);
xcc = xc(:,1);
ycc = yc(1,:);

y_throat = 57;

z = permute(-1*vs_let(trim,'map-sed-series',{0},'DPS',{0,0},'quiet'),[2 3 1]);
h = vs_let(trim,'map-series',{0},'S1',{0,0},'quiet');

rundir=[runsdir filesep runname '_' num2str(idx2,'%1.0f')];
trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
z_2 = permute(-1*vs_let(trim,'map-sed-series',{0},'DPS',{0,0},'quiet'),[2 3 1]);
h_2 = vs_let(trim,'map-series',{0},'S1',{0,0},'quiet');



fig = figure('color','white','Name','Breach_WaterSurface','visible','off','Units','Pixels','PaperUnits','Centimeters','Units','Centimeter',...
    'InvertHardCopy', 'off','NextPlot','add','Position',[5 5 35 15]);

a = tight_subplot(1,2,0.02,0.02);

for ii=1:length(a), hold(a(ii),'on'), box(a(ii),'on'), end

cmap = demcmap([-3 3]);

title(a(1),'Washover')
title(a(2),'Breach')

for kk=1:length(t)
    zmap = squeeze(z(:,:,kk));
    hmap = squeeze(h(kk,:,:));   
    nanmap = (-0.2+hmap)<=zmap;
    hmap = min(-0.1,1+(hmap.*-1));
    hmap(nanmap) = nan;
          
    h0 = pcolor(a(1),xc,yc,zmap);
    h1 = pcolor(a(1),xc,yc,hmap);
    
    shading(a(1),'flat'),
    set(a(1),'CLim',[-3 3]);
    set(a(1),'XLim',[0 500],'YLim',[-200 200],'DataAspectRatio',[1 1 1])
    set(a(1),'Layer','top','XTickLabels','','YTickLabels','')
    
    colormap(cmap)
    cb = colorbar(a(1),'location','east');
    cb.TickLabels = {'3','2','1','0','1','2','3'};
    xlabel(cb,'<- Water elevation (m) | Barrier elevation (m) ->')
    
    
    zmap = squeeze(z_2(:,:,kk));
    hmap = squeeze(h_2(kk,:,:));   
    nanmap = (-0.2+hmap)<=zmap;
    hmap = min(-0.1,1+(hmap.*-1));
    hmap(nanmap) = nan;
          
    h2 = pcolor(a(2),xc,yc,zmap);
    h3 = pcolor(a(2),xc,yc,hmap);
    
    shading(a(2),'flat'),
    set(a(2),'CLim',[-3 3]);
    set(a(2),'XLim',[0 500],'YLim',[-200 200],'DataAspectRatio',[1 1 1])
    set(a(2),'Layer','top','XTickLabels','','YTickLabels','')
   

    
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),64);
    if kk == 1;
        imwrite(imind,cm,savename,'gif', 'Loopcount',inf,'DelayTime',1/20);
    else
        imwrite(imind,cm,savename,'gif','WriteMode','append','DelayTime',1/20);
    end
    
    delete(h0)
    delete(h1)
    delete(h2)
    delete(h3)
    
    
end

imwrite(imind,cm,savename,'gif','WriteMode','append','DelayTime',1/5);