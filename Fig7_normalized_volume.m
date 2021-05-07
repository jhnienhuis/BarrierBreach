clr
model_runs = {'GapSize','Roughness','PeakRoughness','Duration','Width','StormPeak','vegdensity'};
Rsed = 1.65;
gravity = 10;
d50 = 2e-4;
porosity=0.4;
%plot([1 1e9],[1 1e9],'--k'), hold on
%cf=0.1;
Cdrag_bap = @(Cb,Cd,m,d,h,frac) (gravity/Cb^2+frac*Cd*m*d*h/2); %eq. 9.328 Delft Flow Manual, rewritten into a drag coefficient to be: Cdrag = g/Cb^2+Cd*n*h/2, including fraction coverage for the vegetation part

%comparison w/ analytical
Va = @(T,g,s,h,C,R,d50,w,lw) (lw.*(T.*sqrt(2).*sqrt(g).*s.^(5/2).*(s-h).^(5/2))./(560.*C.*R.^2.*d50.*w.^(5/2)));


for ii=1:length(model_runs),
    
    load(model_runs{ii},'output');
    
    if ~isfield(output,'vegdensity'),
        output.vegdensity=ones(size(output.overwash_volume)).*5;
    end
        
    for jj=1:length(output.overwash_volume),
        
        storm_surge = [1 output.storm_peak(jj) output.storm_peak(jj) 1 1]-1;
        storm_time = [0 (12-output.duration(jj))*60 1+(12+output.duration(jj))*60 24*60 26*60].*60;
        Cdrag = Cdrag_bap(45,1,output.vegdensity(jj)./5e-3,5e-3,0.5,output.roughness(jj));
        
        E = predict_volume(storm_time,storm_surge,output.storm_peak(jj),output.height(jj),Cdrag,gravity,Rsed,output.width(jj),d50);
        Vpot{ii}(jj) = E.*output.lipwidth(jj)./2*output.width(jj);
    end
    
    B{ii} = output.breach';
    %output.overwash_volume(output.overwash_volume<1e-4) = 1;
    Vobs{ii} = output.overwash_volume;
    Vbar{ii} = output.vbarrier;
    Vq{ii} = output.qs_cum.*(1/(1-porosity));
    Mb{ii} = output.modelbr;

end
Mb = [Mb{:}]>1;
B = [B{:}];
Vobs = [Vobs{:}];
Vbar = [Vbar{:}];
Vpot = [Vpot{:}];

Vq = cumsum(cat(1,Vq{:}),2);
Vnorm =  Vq'./Vbar;
t = linspace(0,24,size(Vnorm,1)); %hrs

subplot(1,2,1)
plot(t,(Vnorm(:,~(B|Mb))),'b'), hold on,
plot(t,Vnorm(:,B|Mb),'r')
ylim([0 5])
xlabel('time (hr)')
xlim([0 24]), set(gca,'XTick',[0 6 12 18 24])
ylabel('observed storm impact')

subplot(2,2,2)
x = (max(1e-3,(Vpot./Vbar)));
y = Vnorm(end,:);
y(y<1e-4) = nan;
%R2 of trend: corrcoef(log10(x),log10(y)).^2 = 0.81
scatter(x,y), hold on
xlim([1e-3,1e2])
scatter(x(B|Mb),y(B|Mb),'r')
set(gca,'XScale','log','YScale','log')
box('on')
xlabel('potential storm impact')
ylabel('observed storm impact')

subplot(2,2,4)
c = [logspace(-3,1,9)';50];
breach_prob = accumarray(discretize(x,c)',B|Mb,size(c),@mean);
plot(log10(c(1:end-1))+0.25,breach_prob(1:end-1),'-o')

xlim([-3 2]), ylim([0 1])
ylabel('fraction breached')

saveas(gcf,'D:\Dropbox\2021 BarrierBreach JGR\Fig7_normalized_volume.svg')