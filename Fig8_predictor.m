clr
model_runs = {'GapSize','Roughness','PeakRoughness','Duration','Width','StormPeak'};
Rsed = 1.65;
gravity = 10;
d50 = 2e-4;
plot([1 1e9],[1 1e9],'--k'), hold on
cf=0.1;

%comparison w/ analytical
Va = @(T,g,s,h,C,R,d50,w,lw) (lw.*(T.*sqrt(2).*sqrt(g).*s.^(5/2).*(s-h).^(5/2))./(560.*C.*R.^2.*d50.*w.^(5/2)));


for ii=1:length(model_runs),
    
    load(model_runs{ii},'output');
    %V{ii} = zeros(size(output.overwash_volume));
    for jj=1:length(output.overwash_volume),
        
        storm_surge = [1 output.storm_peak(jj) output.storm_peak(jj) 1 1]-1;
        storm_time = [0 (12-output.duration(jj))*60 1+(12+output.duration(jj))*60 24*60 26*60].*60;
        Cdrag{ii}(jj) = gravity./(55-(30*(output.roughness(jj)))).^2;
        E = predict_volume(storm_time,storm_surge,output.storm_peak(jj),output.height(jj),Cdrag{ii}(jj),gravity,Rsed,output.width(jj),d50).*cf;
        Vowt{ii}(jj) = E.*output.lipwidth(jj)./2*output.width(jj);
    end
    
    Vap{ii} = Va(output.duration.*3600,gravity,output.storm_peak,output.height,gravity./(55-(30*(output.roughness))).^2,Rsed,d50,output.width,output.lipwidth);
    B{ii} = output.breach';
    output.overwash_volume(output.overwash_volume<1) = 1;
    Vowd3d{ii} = output.overwash_volume;
    Vbar{ii} = output.vbarrier; %height.*output.lipwidth.*output.width;
    Mb{ii} = output.modelbr;
    Bw{ii} = output.width;
    Gh{ii} = output.height;
    Gw{ii} = output.lipwidth;
    Smax{ii} = output.storm_peak;
    Tstorm{ii} = output.duration;
    Ow_width{ii} = output.overwash_width;
    Ow_length{ii} = output.overwash_length;
    
end
Vowt = [Vowt{:}];
Vowd3d = [Vowd3d{:}];
Vbar = [Vbar{:}];
Mb = [Mb{:}];
B = [B{:}];
B(Mb>1) = 1; %count model breakdowns as breaches
subplot(1,2,1)
scatter(Vowt(~B),Vowd3d(~B),'ob','filled')
set(gca,'XScale','log','YScale','log'), hold on
xlim([1e1 1e6]),ylim([1e-1 1e6]);
%set(gca,'XTick',[logspace(-2, 8,6)])

box('on')
xlabel('Potential overwash volume (m3)')
ylabel('Observed overwash volume (m3)')
%scatter(Vpot(B),Vobs(B),'or','filled')
scatter(Vowt(B),1e5*(ones(sum(B),1)-0.5+rand(sum(B),1)),'ob','filled')

%for supplementary table:
t = table((1:length(B))',[Bw{:}]',[Gh{:}]',[Gw{:}]',[Smax{:}]',[Tstorm{:}]',[Cdrag{:}]',double(B)',[Ow_width{:}]',[Ow_length{:}]',Vowd3d',Vowt',Vbar');

%% add sandy data
sandy = xlsread('Nienhuis_BarrierBreach_Tables.xlsx','sandy');
a = sandy(:,7)==1;

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;
lipwidth=100;
lip_height = max(0.3,sandy(:,5)-(2*sandy(:,6)));
width=sandy(:,4);
h_limiter = 0.3;
Vbar_sandy = lip_height.*lipwidth.*width;
s = size(sandy,1);

%


Vowt1 = zeros(s,1);
for ii=1:s,
storm_time = sandy(ii,9)*3600.*[0 0.5 1 1.5 2];
storm_surge = [0 sandy(ii,8) sandy(ii,8) 0 0];

%Cdrag = gravity./((55-(sandy(ii,7).*40)).^2);
Cdrag = gravity./((65-(sandy(ii,7).*50)).^2);

E = predict_volume(storm_time,storm_surge,sandy(ii,8),min(sandy(ii,8)-0.3,lip_height(ii)),Cdrag,gravity,Rsed,sandy(ii,4),d50);

Vowt1(ii) = E.*lipwidth*sandy(ii,4);
%Eobs = max(z(22,2:end-1,end))
end

%https://cera.coastalrisk.live/s/e2e0
overwash_volume = sandy(:,11)*0.3; %from eli's data
Vowobs = overwash_volume;
Vowobs(Vowobs<1) = 1e5; %breach
scatter(Vowt1(a),Vowobs(a),'sr','filled','markeredgecolor','k')
scatter(Vowt1(~a),Vowobs(~a),'sg','filled','markeredgecolor','k')
set(gca,'XScale','log','YScale','log'),% hold on

legend
%set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 8.1, 8.1]);
set(gca, 'FontSize', 8)
%
subplot(1,2,2)
scatter(Vowt(~B)./Vbar(~B),Vowd3d(~B)./Vbar(~B),'ob','filled'), hold on

scatter(Vowt(B)./Vbar(B),10*(ones(sum(B),1)-0.5+rand(sum(B),1)),'ob','filled')
set(gca,'XScale','log','YScale','log'), hold on
Vnormobs = overwash_volume./Vbar_sandy;
Vnormobs(overwash_volume<1) = 10; %breach

scatter(Vowt1(a)./Vbar_sandy(a),Vnormobs(a),'sr','filled','markeredgecolor','k')
scatter(Vowt1(~a)./Vbar_sandy(~a),Vnormobs(~a),'sg','filled','markeredgecolor','k')

%percentage of breaches for pot>1
pred = Vowt1(~a)./Vbar_sandy(~a);
obs = Vnormobs(~a);
sum(pred>1 & obs>1)./sum(pred>1);
%percentage of overwashes for pot<1
sum(pred<1 & obs<1)./sum(pred<1);


xlim([1e-3 1e2]),ylim([1e-5 1e2]);
%set(gca,'XTick',[logspace(-2, 1,6)])
box('on')
xlabel('Potential overwash volume (m3)')
ylabel('Observed overwash volume (m3)')
set(gca, 'FontSize', 8)
saveas(gcf,'Fig8_predictor.svg')

%supplementary table for sandy
t = table(Vowobs,Vowt1,Vbar_sandy);
