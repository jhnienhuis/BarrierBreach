function run_series_of_models

%we want to predict breach vs. overwash deposit
%run series of simulations
%comparison against sandy

%vary storm duration/max height
%barrier width/roughness

%extract breach yes/no
%volume overwash
%

%which parameters

out.storm_peak = 3; %3 is default
out.width = 300; %300 is default
out.height = 1:0.5:2.5; %lip height, 2 is default
out.roughness = 0.1; %10% is default. percentage coverage  vegetation on island(?) density and height: https://link.springer.com/content/pdf/10.1007/s11258-005-2467-5.pdf
out.duration = 1; %1 is default
out.lipwidth = 10:10:100; %25 is default

%put parameter in struct
out_names = fieldnames(out);
out_cell = struct2cell(out);
sizes = cellfun('prodofsize', out_cell);

%other relevant parameters
out.runsdir = 'F:\Delft3D\BarrierIsland';
out.runname = ['Test'];
runid='bypass';

inputdir=[out.runsdir filesep 'ModelRun1' filesep];



for irun = 1:prod(sizes)
    [a] = ind2subv(sizes,irun);
    [p] = struct_fun(out_names,out_cell,a);
    
    %adjust duration
    storm_surge = [1 p.storm_peak p.storm_peak 1 1];
    storm_time = [0 (12-p.duration)*60 1+(12+p.duration)*60 24*60 26*60];
    
    rundir=[out.runsdir filesep out.runname '_' num2str(irun,'%1.0f')];
    mkdir(rundir);
    copyfile([inputdir '*.*'],rundir);
    
       
    %adjust storm_peak
    for tt=1:5, 
        findreplace([rundir filesep runid '.bct'],['t.0' num2str(tt) '00000'],[num2str(storm_time(tt),'%1.7f')]);
        findreplace([rundir filesep runid '.bct'],[char(tt+96) '.0000000'],[num2str(storm_surge(tt),'%1.7f')]); 
    end
    
    %adjust width
    island_edge = dep_maker_func(rundir,p.width,p.height,p.lipwidth);
    
    %adjust vegetation
    findreplace([rundir filesep runid '.veg'],['1 1 171 111 1 0.8'],['1 1 171 111 1 ' num2str(1-p.roughness,'%1.3f')]); 
    findreplace([rundir filesep runid '.veg'],['10 1 60 111 2 0.2'],['10 1 ' num2str(island_edge,'%2.0f') ' 111 2 ' num2str(p.roughness,'%1.3f')]); 

    
    
    % And run the simulation
    curdir=pwd;
    cd(rundir);
    dos('call run_flow2d3d_parallel.bat'); % must sit in inputdir!
    cd(curdir);
    
end

save([out.runname '.mat'],'-struct','out')


end

function [p] = struct_fun(q_names,q_cell,a)

for j=1:length(q_names), p.(q_names{j}) = q_cell{j}(a(j));  end

end

function v = ind2subv(siz,ndx)
[out{1:length(siz)}] = ind2sub(siz,ndx);
v = cell2mat(out);
end










