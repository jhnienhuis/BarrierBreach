clr
model_runs = {'GapSize','Roughness','PeakRoughness','Duration','Width','StormPeak','vegdensity'};
Rsed = 1.65;
gravity = 10;
d50 = 2e-4;
plot([1 1e9],[1 1e9],'--k'), hold on
Cdrag_bap = @(Cb,Cd,m,d,h,frac) (gravity/Cb^2+frac*Cd*m*d*h/2); %eq. 9.328 Delft Flow Manual, rewritten into a drag coefficient to be: Cdrag = g/Cb^2+Cd*n*h/2, including fraction coverage for the vegetation part

%comparison w/ analytical
Va = @(T,g,s,h,C,R,d50,w,lw) (lw.*(T.*sqrt(2).*sqrt(g).*s.^(5/2).*(s-h).^(5/2))./(560.*C.*R.^2.*d50.*w.^(5/2)));


for ii=1:length(model_runs),
    
    load(model_runs{ii},'output');
    %if ~isfield(output,'vegdensity'),
    %    output.vegdensity=ones(size(output.overwash_volume))*5;
    %end
    %V{ii} = zeros(size(output.overwash_volume));
    for jj=1:length(output.overwash_volume),
        
        storm_surge = [1 output.storm_peak(jj) output.storm_peak(jj) 1 1]-1;
        storm_time = [0 (12-output.duration(jj))*60 1+(12+output.duration(jj))*60 24*60 26*60].*60;
        Cdrag{ii}(jj) = Cdrag_bap(45,1.65,output.vegdensity(jj),0.5,0.5,output.roughness(jj));
        E = predict_volume(storm_time,storm_surge,output.storm_peak(jj),output.height(jj),Cdrag{ii}(jj),gravity,Rsed,output.width(jj),d50);
        Vowt{ii}(jj) = E.*output.lipwidth(jj)./2*output.width(jj);
    end
    
    Vap{ii} = Va(output.duration.*3600,gravity,output.storm_peak,output.height,Cdrag_bap(45,1,output.vegdensity(jj),1,0.5,output.roughness(jj)),Rsed,d50,output.width,output.lipwidth);
    B{ii} = output.breach';
    output.overwash_volume(output.overwash_volume<1) = 1e-1;
    Vowd3d{ii} = output.overwash_volume;
    Vbar{ii} = output.vbarrier;
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
Vowd3d(Vowd3d<1) = nan;
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
xlabel('Predicted overwash volume (m3)')
ylabel('Observed overwash volume (m3)')
%scatter(Vpot(B),Vobs(B),'or','filled')
scatter(Vowt(B),1e5*(ones(sum(B),1)-0.5+rand(sum(B),1)),'ob','filled')

%for supplementary table:
t = table((1:length(B))',[Bw{:}]',[Gh{:}]',[Gw{:}]',[Smax{:}]',[Tstorm{:}]',[Cdrag{:}]',double(B)',[Ow_width{:}]',[Ow_length{:}]',Vowd3d',Vowt',Vbar');

%% add sandy data
sandy = xlsread('Nienhuis_BarrierBreach_Tables.xlsx','sandy');
a = sandy(:,7)>1e-1;

d50 = 2e-4;
gravity = 10;
Rsed = 1.65;
lipwidth=50;
h_limiter = 0.3;

lip_height = max(h_limiter,sandy(:,5)-(2*sandy(:,6)));
width=sandy(:,4);
Vbar_sandy = lip_height.*lipwidth.*width;
s = size(sandy,1);

Vowt1 = zeros(s,1);

for ii=1:s,
storm_time = sandy(ii,9)*3600.*[0 0.5 1 1.5 2];
storm_surge = [0 sandy(ii,8) sandy(ii,8) 0 0];

E = predict_volume(storm_time,storm_surge,sandy(ii,8),min(sandy(ii,8)-h_limiter,lip_height(ii)),sandy(ii,7),gravity,Rsed,sandy(ii,4),d50);

Vowt1(ii) = E.*lipwidth*sandy(ii,4);

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
xlabel('Predicted overwash volume (m3)')
ylabel('Observed overwash volume (m3)')
set(gca, 'FontSize', 8)
saveas(gcf,'D:\Dropbox\2021 BarrierBreach JGR\Fig8_predictor.svg')

%supplementary table for sandy
t = table(Vowobs,Vowt1,Vbar_sandy);
