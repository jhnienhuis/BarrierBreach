clr
Rsed = 1.65;
gravity = 10;
d50 = 2e-4;
porosity=0.4;
%plot([1 1e9],[1 1e9],'--k'), hold on
cf=0.1;
    
    load('Width','output');

    for jj=1:length(output.overwash_volume),
        
        storm_surge = [1 output.storm_peak(jj) output.storm_peak(jj) 1 1]-1;
        storm_time = [0 (12-output.duration(jj))*60 1+(12+output.duration(jj))*60 24*60 26*60].*60;
        Cdrag = gravity./(55-(30*(output.roughness(jj)))).^2;
        E = predict_volume(storm_time,storm_surge,output.storm_peak(jj),output.height(jj),Cdrag,gravity,Rsed,output.width(jj),d50).*cf;
        Vpot(jj) = E.*output.lipwidth(jj)./2*output.width(jj);
    end
    
    
    B = output.breach';
    
    
    
    Vobs = output.overwash_volume;
    Vbar = output.height.*output.lipwidth.*output.width./2;
    Vq = sum(cat(1,output.qs_cum.*(1/(1-porosity))),2);
    Mb = output.modelbr;
    
subplot(1,2,1)
plot(output.width,output.overwash_length,'-o'), hold on,
plot(output.width+output.overwash_length,output.overwash_length,'-o')
xlabel('Barrier Width (m)')
ylabel('Washover distance (m)')
set(gca, 'FontSize', 7,'FontName','Helvetica')
subplot(1,2,2)
plot(output.width,Vq./output.lipwidth,'-o'), hold on
plot(output.width,output.width.*output.height,'-o')
xlabel('Barrier Width (m)')
ylabel('Alongshore average Overwash flux (m3/m)')

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 7,'FontName','Helvetica')
saveas(gcf,'Fig9_critical_width.svg')