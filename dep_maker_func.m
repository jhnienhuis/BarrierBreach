function [island_edge] = dep_maker_func(direc,width,height,lipwidth)

% profile units in meters
x_profile = 0:2000;

z_profile = -[linspace(0,4,100) linspace(4,0,width-100) linspace(0,-3,300) -3*ones(1,1701-width)]; %depth is positive in d3d

%% flow
[x,y] = meshgrid([0:5:395 400:10:600 620:20:2000],[500:-20:220,200:-5:-200,-220:-20:-500]);
name = [direc filesep 'grid_flow'];

%dep stuff
C{1} = {'Coordinate';'System';'=';'Cartesian';'999';'999';'0';'0';'0'};
[m,n] = size(y);
C{1}{5} = num2str(m);
C{1}{6} = num2str(n);

%dep file

z = interp1(x_profile,z_profile,x(10,:));

z = repmat(z,m,1);

z(m+1,:) = -999;
z(:,n+1) = -999;

island_edge = find(z(2,:)>0,1)-1;
%% write dep

clear i
fid = fopen(name,'w+');

for i=1:numel(z(1,:)),
    fprintf(fid,'%f ',z(:,i));
    fprintf(fid,'\r\n');
end

fclose(fid)

movefile(name,[name '.dep'])
pause(0.1)

%% add channel to dep file

fid = fopen([name '.dep']);

C = textscan(fid, '%f');
z = reshape(C{1},m+1,n+1);

fclose(fid);

breach_elevation = -1*height; %m, depth is still positive in d3d, so elevation is negative
breach_width = lipwidth; %m
breach_width_cell = round(breach_width./5);
breach_mid_idx = round(size(z,1)./2 + (0:breach_width_cell) - 0.5*breach_width_cell);

z_new = z;
z_new(breach_mid_idx,:) = max(z_new(breach_mid_idx,:),breach_elevation);
z(m+1,:) = -999;
z(:,n+1) = -999;
%pcolor(x(1,:),y(:,1),z_new(1:end-1,1:end-1))

%z_new([1 m+1],:) = -5; add a wall to either side

fid = fopen([name '_breach'],'w+');

for i=1:numel(z_new(1,:)),
    fprintf(fid,'%f ',z_new(:,i));
    fprintf(fid,'\r\n');
end

fclose(fid);

movefile([name '_breach'],[name 'breach.dep']);
drawnow