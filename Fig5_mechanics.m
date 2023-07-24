clr
addpath('D:\Drive\_Tools\delft3d_matlab')
load('PeakRoughness');

%inputdir=[runsdir filesep 'ModelRun1' filesep];

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

%idx = 4; %splay
idx = 27; %breach



%adjust duration
storm_surge = [1 output.storm_peak(idx) output.storm_peak(idx) 1 1];
storm_time = [0 (12-output.duration(idx))*60 1+(12+output.duration(idx))*60 24*60 26*60];

rundir=[runsdir filesep runname '_' num2str(idx,'%1.0f')];

%list of options in delft3d


%rundir='PeakRoughness_4';

trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
trih = vs_use([rundir filesep 'trih-bypass.dat'],[rundir filesep 'trih-bypass.def'],'quiet');
t = vs_let(trim,'map-infsed-serie',{0},'MORFT','quiet');

t = vs_let(trim,'map-infsed-serie',{0});

xc=vs_get(trim,'map-const','XCOR','quiet!');
yc=vs_get(trim,'map-const','YCOR','quiet!');

cell_area = cellarea(xc,yc);
xcc = xc(:,1);
ycc = yc(1,:);

y_mid = 57;
x_mid = 1:length(xcc);

zmap = squeeze(-1*vs_let(trim,'map-sed-series',{1},'DPS',{0,0},'quiet'));

x_throat = find(zmap(:,2)>0);
y_throat = find(abs(zmap(20,:)-output.height(idx))<0.1);
y_throat_ext = y_throat(1)+[-10:(10+length(y_throat))];
y_diff = reshape(abs(diff(ycc([y_throat_ext(1)-1 y_throat_ext]))),1, 1, []);


v_total = sum(cell_area(x_throat,y_throat).*zmap(x_throat,y_throat),'all');
output.height(idx).*output.lipwidth(idx).*output.width(idx).*0.5;



z = -1*vs_let(trim,'map-sed-series',{0},'DPS',{x_mid,y_mid},'quiet')';
u = vs_let(trim,'map-series',{0},'V1',{x_mid,y_mid,0},'quiet');
taub = vs_let(trim,'map-series',{0},'TAUETA',{x_mid,y_mid},'quiet'); %N/m2
h = vs_let(trim,'map-series',{0},'S1',{x_mid,y_mid},'quiet');
qs_bed = vs_let(trim,'map-sed-series',{0},'SBVV',{x_mid,y_throat_ext,1},'quiet'); %m2/s
qs_sus = vs_let(trim,'map-sed-series',{0},'SSVV',{x_mid,y_throat_ext,1},'quiet');

tt=round([0.01 0.25 0.5 0.75 1].*length(t));

for ii=1:length(tt),
%x plots
subplot(3,5,ii);
plot(xcc(x_mid),h(tt(ii),:)), hold on,
%plot(xcc(x_throat),u(tt(ii),:)), 
plot(xcc(x_mid),z(:,tt(ii)));
plot(xcc(x_mid),zeros(size(x_mid)),':')
ylim([-3 3])
xlim([0 500])
set(gca,'XTick',[0 250 500])
subplot(3,5,ii+5)
plot(xcc(x_mid),taub(tt(ii),:)), hold on,
plot(xcc(x_mid),qs_bed(tt(ii),:,y_throat_ext==y_mid).*1000), hold on,
plot(xcc(x_mid),qs_sus(tt(ii),:,y_throat_ext==y_mid).*1000), hold on,
ylim([0 100])
xlim([0 500])
set(gca,'XTick',[0 250 500])
xlabel('across (m)')
end





subplot(3,1,3)

%calculate total volume of sediment
qw = vs_let(trih,'his-series',{0},'CTR',{1},'quiet'); %m3s-1
wlhead = (h(:,1)-h(:,end)); %./range(xcc(x_throat)

qs = sum((qs_bed+qs_sus).*y_diff,3);
qs_cum = cumsum(qs(2:end,53)).*diff(t)*3600*24;
plot(t(2:end).*24,qs_cum./v_total), hold on
plot(t.*24,u(:,53))
plot(t.*24,wlhead)
plot(storm_time./60,(storm_surge-1)),
xlabel('time (hr)')
xlim([0 24]), set(gca,'XTick',[0 6 12 18 24])
ylim([0 3.5])


saveas(gcf,'D:\Dropbox\2021 BarrierBreach JGR\Fig5_mechanics.svg')