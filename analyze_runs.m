function analyze_runs

model_runs = {'GapSize','Roughness','PeakRoughness','Duration','Width','StormPeak','vegdensity'};

for ii=1:length(model_runs),
    analyze_run(model_runs{ii});
end
end

function analyze_run(name)

addpath('D:\Dropbox\_Tools\delft3d_matlab')
out = load(name);
if ~isfield(out,'lipwidth'), out.lipwidth=50; end
if ~isfield(out,'vegdensity'), out.vegdensity=5; end
runsdir = out.runsdir; out = rmfield(out,'runsdir');
runname =out.runname; out = rmfield(out,'runname');
if isfield(out,'output'), out = rmfield(out,'output'); end
inputdir=[runsdir filesep 'ModelRun1' filesep];



%put parameter in struct
out_names = fieldnames(out);
out_cell = struct2cell(out);
sizes = cellfun('prodofsize', out_cell);

%other relevant parameters


d50 = 2e-4;
gravity = 10;
Rsed = 1.65;

z = zeros(172,112,prod(sizes));

for irun = 1:prod(sizes)
    irun
    [a] = ind2subv(sizes,irun);
    [p] = struct_fun(out_names,out_cell,a);
    
    %adjust duration
    storm_surge = [1 p.storm_peak p.storm_peak 1 1];
    storm_time = [0 (12-p.duration)*60 1+(12+p.duration)*60 24*60 26*60];
    
    rundir=[runsdir filesep runname '_' num2str(irun,'%1.0f')];
    
    
  
    trim = vs_use([rundir filesep 'trim-bypass.dat'],[rundir filesep 'trim-bypass.def'],'quiet');
    t = vs_let(trim,'map-infsed-serie',{0},'MORFT','quiet');
    
    %roughness_map = vs_let(trim,'map-series',{0},'CFUROU',{0,0},'quiet');
    
    if irun==1,
    xc=vs_get(trim,'map-const','XCOR','quiet!');
    yc=vs_get(trim,'map-const','YCOR','quiet!');
    ycc = yc(1,:);

    cell_area = cellarea(xc,yc);
    end
    
    z(:,:,irun) = squeeze(vs_let(trim,'map-sed-series',{length(t)},'DPS',{0,0},'quiet'));
    z0 = squeeze(vs_let(trim,'map-sed-series',{1},'DPS',{0,0},'quiet'));
        
    idx = find(z0(:,2)>0,1):170;
    idy = 2:111; %57+[-8:8];
    
    dz = z0-z(:,:,irun);

    dz(abs(dz)>10) = 0;
    
    
    x_throat = find(z0(:,2)<0,1,'last');
    y_throat = find(abs(z0(20,:)+p.height)<0.1);
    y_throat_ext = y_throat(1)+[-10:(10+length(y_throat))];
    y_diff = reshape(abs(diff(ycc([y_throat_ext(1)-1 y_throat_ext]))),1, 1, []);
    
    qs_bed = vs_let(trim,'map-sed-series',{0},'SBVV',{x_throat,y_throat_ext,1},'quiet');
    qs_sus = vs_let(trim,'map-sed-series',{0},'SSVV',{x_throat,y_throat_ext,1},'quiet');
    qs = [0;sum((qs_bed(2:end,:,:)+qs_sus(2:end,:,:)).*diff(t)*3600*24.*y_diff,[3])]; %m3

    output.modelbr(irun) = max(abs(diff(dz)),[],'all');
    output.vbarrier(irun) = sum(cell_area(x_throat,y_throat).*z0(1:x_throat,y_throat),'all').*-1;
    output.overwash_volume(irun) = squeeze(sum(cell_area(idx,idy).*dz(idx,idy),[1 2]));
    output.qs_cum(irun,1:length(t)) = qs;
    output.width(irun) = p.width;
    output.storm_peak(irun) = p.storm_peak;
    output.roughness(irun) = p.roughness;
    output.height(irun) = p.height;
    output.lipwidth(irun) = p.lipwidth;
    output.duration(irun) = p.duration;
    output.vegdensity(irun) = p.vegdensity;
    
    [col,row] = find(z(idx,idy,irun)<0);
    if isempty(col),
        output.overwash_length(irun) = 0;
        output.overwash_width(irun) = 0;
    else,
        output.overwash_length(irun) = xc(max(col)+idx(1)-1)-xc(idx(1),1);
        output.overwash_width(irun) = yc(1,min(row)+idy(1)-1)-yc(1,max(row)+idy(1)-1);
    end
    
    %figure,
    %imagesc(z(:,:,irun),[-3 3])
    %grid on
    %pause(0.3)
    %squeeze(max(z(22,:,:),[],2));

end
  
output.breach = squeeze(max(z(22,:,:),[],2))>0.01; %squeeze(mean(z(20:40,57,:)))>0; %

output.overwash_volume(output.overwash_volume<1e-3) = 1e-3;

%{
for ii=[8 16 20],
    figure,
    imagesc(z(:,:,ii),[-3 3])
    grid on
    pause(0.3)
squeeze(max(z(22,:,:),[],2));

end
%}

save(name,'output','-append')

end


function [p] = struct_fun(q_names,q_cell,a)

for j=1:length(q_names), p.(q_names{j}) = q_cell{j}(a(j));  end

end

function v = ind2subv(siz,ndx)
[out{1:length(siz)}] = ind2sub(siz,ndx);
v = cell2mat(out);
end
