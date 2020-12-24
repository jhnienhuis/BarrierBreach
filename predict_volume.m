function E = predict_volume(storm_time,storm_surge,storm_peak,height,Cdrag,gravity,Rsed,width,d50)
dz = @(t) (interp1(storm_time,storm_surge,t,'linear','extrap'));
waterdepth = max(0.01,(storm_peak-height))./2;
Ez = @(t) (0.05/Cdrag.*waterdepth.^2.5*dz(t).^2.5.*sqrt(gravity)./d50./Rsed.^2./width.^3.5);
E = integral(Ez,0,1e5);
end
