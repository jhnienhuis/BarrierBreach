clr
model_runs = {'Roughness','Duration','Width','StormPeak','GapSize','GapSize','vegdensity'};
mid_par = [6,6,7,11,3,17,6];
par = {'roughness','duration','width','storm_peak','height','lipwidth','vegdensity'};

for ii=1:length(model_runs),
 

    load(model_runs{ii});
    
    if ii==5, 
        i1 = 1; i2 = 4;
    elseif ii==6, 
        i1 = 4; i2 = 40;
    else
        i1=1; i2=length(eval(par{ii}));
    end
    
    x = output.(par{ii})(1:i1:i2)./output.(par{ii})(mid_par(ii));
    
    
    plot(x,output.overwash_volume(1:i1:i2)./output.overwash_volume(mid_par(ii)),'-o')
    
    hold on, grid on    

    
    
end

set(gca,'YScale','log')

ylabel('Overwash volume (compared to default estimate)')
xlabel('Parameter (compared to default estimate)')
legend('Vegetation cover','Storm duration','Barrier width','Storm surge height','Dune gap elevation','Dune gap width','Vegetation density')

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'D:\Dropbox\2021 BarrierBreach JGR\barrier_breach_data\FigS1_sensitivity.svg')

%figure shows the relative change of overwashing volumes compared to the
%relative change of various Delft3D input parameters. Overwash volumes are 
%most sensitive to storm surge height, barrier width, and barrier elevation
%(dune gap) are most