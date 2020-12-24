clr
model_runs = {'GapSize','Roughness','PeakRoughness','Duration','Width','StormPeak'};
Rsed = 1.65;
gravity = 10;
d50 = 2e-4;
porosity=0.4;
%plot([1 1e9],[1 1e9],'--k'), hold on
cf=0.1;

%comparison w/ analytical
Va = @(T,g,s,h,C,R,d50,w,lw) (lw.*(T.*sqrt(2).*sqrt(g).*s.^(5/2).*(s-h).^(5/2))./(560.*C.*R.^2.*d50.*w.^(5/2)));


for ii=1:length(model_runs),
    
    load(model_runs{ii},'output');
    %V{ii} = zeros(size(output.overwash_volume));
    for jj=1:length(output.overwash_volume),
        
        storm_surge = [1 output.storm_peak(jj) output.storm_peak(jj) 1 1]-1;
        storm_time = [0 (12-output.duration(jj))*60 1+(12+output.duration(jj))*60 24*60 26*60].*60;
        Cdrag = gravity./(55-(30*(output.roughness(jj)))).^2;
        E = predict_volume(storm_time,storm_surge,output.storm_peak(jj),output.height(jj),Cdrag,gravity,Rsed,output.width(jj),d50).*cf;
        Vpot{ii}(jj) = E.*output.lipwidth(jj)./2*output.width(jj);
    end
    
    Vap{ii} = Va(output.duration.*3600,gravity,output.storm_peak,output.height,gravity./(55-(30*(output.roughness))).^2,Rsed,d50,output.width,output.lipwidth);
    B{ii} = output.breach';
    output.overwash_volume(output.overwash_volume<1e-3) = 1e-3;
    Vobs{ii} = output.overwash_volume;
    %Vbar{ii} = output.height.*output.lipwidth.*output.width;
    Vbar{ii} = output.vbarrier;
    Vq{ii} = output.qs_cum.*(1/(1-porosity));
    Mb{ii} = output.modelbr;

end
Mb = [Mb{:}]>1;
B = [B{:}];
Vobs = [Vobs{:}];
Vbar = [Vbar{:}];
Vap = [Vap{:}];
Vpot = [Vpot{:}];

Vq = cumsum(cat(1,Vq{:}),2);
Vnorm =  Vq'./Vbar;
t = linspace(0,24,length(Vnorm)); %hrs

subplot(1,2,1)
plot(t,(Vnorm(:,~(B|Mb))),'b'), hold on,
plot(t,Vnorm(:,B|Mb),'r')
ylim([0 5])
xlabel('time (hr)')
xlim([0 24]), set(gca,'XTick',[0 6 12 18 24])
ylabel('observed storm impact')

subplot(2,2,2)
x = (max(1e-2,(Vpot./Vbar)));
y = (max(1e-4,Vnorm(end,:)));
scatter(x,y), hold on
%xlim([1e-2,1e1])
scatter(x(B|Mb),y(B|Mb),'r')
set(gca,'XScale','log','YScale','log')
box('on')
xlabel('potential storm impact')
ylabel('observed storm impact')

subplot(2,2,4)
c = logspace(-2,1,7)';
breach_prob = accumarray(discretize(x,c)',B|Mb,size(c),@mean);
%bar(log10(c)+0.25,[breach_prob 1-breach_prob],'stacked'),
plot(log10(c(1:end-1))+0.25,breach_prob(1:end-1),'-o')

%Vbreach = c(breach_prob>0.5);
xlim([-2 1]), ylim([0 1])
ylabel('fraction breached')
%plot idealized:
%phase space of storm water slope

saveas(gcf,'Fig7_normalized_volume.svg')