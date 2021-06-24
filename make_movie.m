%% breach movie
addpath('D:\Dropbox\_Tools\delft3d_matlab')
load('PeakRoughness');

inputdir=[runsdir filesep 'ModelRun1' filesep];

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

idx = 27; %breach
%idx = 4; %overwash


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

%calculate total volume of sediment
%qw = vs_let(trih,'his-series',{0},'CTR',{1},'quiet'); %m3s-1
%qs_bed = vs_let(trih,'his-sed-series',{0},'SBTR',{1,1},'quiet');
%qs_sus = vs_let(trih,'his-sed-series',{0},'SSTR',{1,1},'quiet'); %kgm-3

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








%{
fig = figure('color','white','Name','Vegetation+Compaction','visible','off','Units','Pixels','PaperUnits','Centimeters','Units','Centimeter',...
    'InvertHardCopy', 'off','NextPlot','add','Position',[5 5 45 15]);

a = tight_subplot(1,2,0.1,0.1,0.1);

tick_max = max([splay_v./1e6;s_discharge./1e3]);
dis_tick = (10^floor(log10(tick_max)))*ceil(tick_max/(10^floor(log10(tick_max))));


set(a(1),'CLim',[-2 4],'NextPlot','add','XLim',[-0.25 10],'YLim',[-3 3],'Layer','top','DataAspectRatio',[1 1 1])
set(a(2),'NextPlot','add','XLim',[0 10],'YLim',[0 dis_tick],'Layer','top')
set(a(2),'PlotBoxAspectRatio',get(a(1),'PlotBoxAspectRatio'))

c1 = colorbar('peer',a(1),'East','color','w');
xlabel(c1,'Bed level (m)','color','w')
box(a(2),'on')
box(a(1),'on')

xlabel(a(1),'Cross-levee (km)');
ylabel(a(1),'Along-levee (km)');

xlabel(a(2),'Time (yr)');
ylabel(a(2),{'Discharge into Splay (10^3 m^3s^-^1)';'\color{red}Splay Volume (10^6 m^3)'});

xticklabels(a(1),'auto')
xticklabels(a(2),'auto')
yticklabels(a(1),'auto')
yticklabels(a(2),'auto')

for i=1:numel(t)
    
    i_dis = ceil(i*length(t_discharge)./length(t));
    
    if mod(i,10)==1,i, end
    
    title(a(1),[' Time: ' num2str(t(i)./365,'%3.1f') ' years'],'interpreter','none')
    
    %h2 = plot(a(1),x(2:end),[-1*zchannel(i,(2:end)); wlchannel(i,(2:end)); u(i,2:end)]);
    
    
    h1 = pcolor(x./1000,y(2:end)./1000,squeeze(-1*bl(i,:,2:end))','Parent',a(1)); shading(a(1),'flat')
    
    h2 = plot(a(2),t_discharge./365, s_discharge./1e3,'-k',t_discharge(i_dis)./365, s_discharge(i_dis)./1e3,'ok');
    h3 = plot(a(2),t_discharge./365, splay_v./1e6,'-r',t_discharge(i_dis)./365, splay_v(i_dis)./1e6,'or');
    
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),64);
    if i == 1;
        imwrite(imind,cm,[case_name '.gif'],'gif', 'Loopcount',inf,'DelayTime',1/30);
    else
        imwrite(imind,cm,[case_name '.gif'],'gif','WriteMode','append','DelayTime',1/30);
    end
    
    delete(h1)
    delete(h2)
    delete(h3)
    
end

imwrite(imind,cm,[case_name '.gif'],'gif','WriteMode','append','DelayTime',1/5);


%% closure movie
case_name = 'BatchMidDomain_Head5_trial12';
case_name2 = 'BatchMidDomain_Veg3';
load([dropbox filesep 'Splay' filesep case_name])
c2 = load([dropbox filesep 'Splay' filesep case_name2]);
y(1) = 3200; y(2) = 3100; y(3) = 3000;

y(end-1) = -3100; y(end) = -3200;
fig = figure('color','white','Name','Vegetation+Compaction','visible','off','Units','Pixels','PaperUnits','Centimeters','Units','Centimeter',...
    'InvertHardCopy', 'off','NextPlot','add','Position',[5 5 15 13.5]);

a = tight_subplot(1,2,0.03,0.03,0.01);
%demcmap([-2.5 3.5])
set(a(1),'CLim',[-2 4],'NextPlot','add','YLim',[-0.25 10],'XLim',[-3 3],'Layer','top','DataAspectRatio',[1 1 1])
set(a(2),'CLim',[-2 4],'NextPlot','add','YLim',[-0.25 10],'XLim',[-3 3],'Layer','top','DataAspectRatio',[1 1 1])
set(a(2),'PlotBoxAspectRatio',get(a(1),'PlotBoxAspectRatio'))

c1 = colorbar('peer',a(2),'North','color','w');
xlabel(c1,'Bed level (m)','color','w')
box(a(2),'on')
box(a(1),'on')

%xlabel(a(1),'Avulsion');
%ylabel(a(1),'Along-levee (km)');
%xlabel(a(2),'Successful diversion');
%ylabel(a(2),'Along-levee (km)');

set(a(1),'fontsize',14)
set(a(2),'fontsize',14)

xticklabels(a(1),'')
xticklabels(a(2),'')
yticklabels(a(1),'')
yticklabels(a(2),'')
%set(gcf,'visible','on')
%
for i=1:3:700
    
    
    
    if mod(i,10)==1,i, end
    
    title(a(1),[' Time: ' num2str(t(ceil(i/10))./365,'%3.1f') ' years'],'interpreter','none')
    
    h1 = pcolor(y(2:end)./1000,x./1000,squeeze(-1*bl(ceil(i/10),:,2:end)),'Parent',a(1)); shading(a(1),'flat')
    %text(a(1),-2.5,10,'A','fontsize',25)
    title(a(2),[' Time: ' num2str(c2.t(i)./365,'%3.1f') ' years'],'interpreter','none')
    h2 = pcolor(y(2:end)./1000,x./1000,squeeze(-1*c2.bl(i,:,2:end)),'Parent',a(2)); shading(a(2),'flat')
    %text(a(2),-2.5,10,'B','fontsize',25)
    %set(gcf,'Visible','on'), break
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),32);
    if i == 1;
        imwrite(imind,cm,[case_name 'comparison2.gif'],'gif', 'Loopcount',inf,'DelayTime',1/30);
    else
        imwrite(imind,cm,[case_name 'comparison2.gif'],'gif','WriteMode','append','DelayTime',1/30);
    end
    
    delete(h1)
    delete(h2)
    
end

imwrite(imind,cm,[case_name 'comparison3.gif'],'gif','WriteMode','append','DelayTime',1/5);
%}