close all
clear all
clc


%% Inputs and outputs

% Input 1 for this script is the bin configuration file which will have a
% name like 'ATGR1_AAF_GV135N33GSIFM.txt'

input_file = 'ATGR1_AAF_GV135N33GSIFM.txt';
output_tx_file_name = 'TxFile.txt';
output_rx_file_name = 'RxFile.txt';
output_tx_to_rx_file_name = 'TxToRxFile.txt';
nTx = 24;

if ~exist(input_file, 'file')
    disp('Could not read file');
    return;
elseif exist(output_tx_file_name, 'file')
    disp([output_tx_file_name ' already exists']);
    return;
elseif exist(output_rx_file_name, 'file')
    disp([output_rx_file_name ' already exists']);
    return;
elseif exist(output_tx_to_rx_file_name, 'file')
    disp([output_tx_to_rx_file_name ' already exists']);
    return;
end

% Input 2 for this script is the direction of current inside the antenna.
% We compute the polarization of the antenna based on this and
% the right hand rule. If we have cables which are mounted on the
% roof and come down the bin wall we expect the polarization to be clockwise.
% Alternatively if the cable it mounted at the floor and goes up towards
% the antenna we expect the polarization to be counter-clockwise
% In receiver mode we always have the opposite orientation of the transmitter

Tx_J_vector = [0 0 -1]; % cable_direction = 'roof_to_floor';
Rx_J_vector = [0 0 +1]; % cable_direction = 'roof_to_floor';

%Tx_J_vector = [0 0 +1]; % cable_direction = 'floor_to_roof';
%Rx_J_vector = [0 0 -1]; % cable_direction = 'floor_to_roof';


%% Read the file and parse data

bin = struct;
antennas = struct;
antenna_start_line = 15; % should be constant

antennas.ideal.Z_mm = [];
antennas.ideal.Phi_rad = [];
antennas.ideal.Coax_mm = [];

antennas.install.Z_mm = [];
antennas.install.Phi_rad = [];
antennas.install.Coax_mm = [];

antennas.raster.Z_mm = [];
antennas.raster.Phi_rad = [];
antennas.raster.Coax_mm = [];

antennas.best.Z_mm = [];
antennas.best.Phi_rad = [];
antennas.best.Coax_mm = [];

fid = fopen(input_file);
disp('Received the following file:');
line = 1;
tline = fgetl(fid);
while ischar(tline)
    disp(['Line #' num2str(line) ': ' tline]);
    tline = fgetl(fid);
    line = line+1;
    
    
    if (line >= 9 && line<=12) % lines with the bin parameters
        parameters = strsplit(tline);
        if line==9
            type = 'ideal'; % ideal positions line
        elseif line==10
            type = 'install'; % install positions line
        elseif line==11
            type = 'raster'; % raster positions line
        elseif line==12
            type = 'best'; % best positions line
        end
        eval(['bin.' type '.height_ft=str2double(parameters{1})']);
        eval(['bin.' type '.plenum_in=str2double(parameters{2})']);
        eval(['bin.' type '.diameter_ft=str2double(parameters{3})']);
        eval(['bin.' type '.angle_deg=str2double(parameters{4})']);
        eval(['bin.' type '.cap_diameter_ft=str2double(parameters{5})']);
        eval(['bin.' type '.num_antenna=str2double(parameters{6})']);
    end
        
    if (line >=antenna_start_line && line <(antenna_start_line+nTx))
        parameters_str = strsplit(tline);
        
        parameters = zeros(13,1);
        for i=1:numel(parameters)
            parameters(i) = str2double(parameters_str{i});
        end
        
        antennas.ideal.Z_mm = [antennas.ideal.Z_mm; parameters(2)];
        antennas.ideal.Phi_rad = [antennas.ideal.Phi_rad; parameters(3)];
        antennas.ideal.Coax_mm = [antennas.ideal.Coax_mm; parameters(4)];
        
        antennas.install.Z_mm = [antennas.install.Z_mm; parameters(5)];
        antennas.install.Phi_rad = [antennas.install.Phi_rad; parameters(6)];
        antennas.install.Coax_mm = [antennas.install.Coax_mm; parameters(7)];
        
        antennas.raster.Z_mm = [antennas.raster.Z_mm; parameters(8)];
        antennas.raster.Phi_rad = [antennas.raster.Phi_rad; parameters(9)];
        antennas.raster.Coax_mm = [antennas.raster.Coax_mm; parameters(10)];
        
        antennas.best.Z_mm = [antennas.best.Z_mm; parameters(11)];
        antennas.best.Phi_rad = [antennas.best.Phi_rad; parameters(12)];
        antennas.best.Coax_mm = [antennas.best.Coax_mm; parameters(13)];
        
    end
end
fclose(fid);


%% Sort through parameters and pull out the ones we want. Conver to SI units
mm_to_m = 0.001;
ft_to_m = 0.3048;
in_to_m = 0.0254;
deg_to_rad = pi/180;

% The prioritization of parameters is best, raster, install, ideal

output.Z_m = NaN*zeros(nTx,1);
idxs = find(isnan(output.Z_m)); % find idxs that are nan
output.Z_m(idxs) = antennas.best.Z_mm(idxs)*mm_to_m; % set output accordingly

idxs = find(isnan(output.Z_m));
output.Z_m(idxs) = antennas.raster.Z_mm(idxs)*mm_to_m; % idxs that are NaN change to raster

idxs = find(isnan(output.Z_m));
output.Z_m(idxs) = antennas.install.Z_mm(idxs)*mm_to_m; % idxs that are NaN change to install

idxs = find(isnan(output.Z_m));
output.Z_m(idxs) = antennas.ideal.Z_mm(idxs)*mm_to_m; % idxs that are NaN change to ideal

% Now do the same for Phi
output.Phi_rad = NaN*zeros(nTx,1);
idxs = find(isnan(output.Phi_rad)); % find idxs that are nan
output.Phi_rad(idxs) = antennas.best.Phi_rad(idxs); % set output to best

idxs = find(isnan(output.Phi_rad));
output.Phi_rad(idxs) = antennas.raster.Phi_rad(idxs); % idxs that are NaN change to raster

idxs = find(isnan(output.Phi_rad));
output.Phi_rad(idxs) = antennas.install.Phi_rad(idxs); % idxs that are NaN change to install

idxs = find(isnan(output.Phi_rad));
output.Phi_rad(idxs) = antennas.ideal.Phi_rad(idxs); % idxs that are NaN change to ideal


% Now do the same for Coax
output.Coax_m = NaN*zeros(nTx,1);
idxs = find(isnan(output.Coax_m)); % find idxs that are nan
output.Coax_m(idxs) = antennas.best.Coax_mm(idxs)*mm_to_m; % set output to best

idxs = find(isnan(output.Coax_m));
output.Coax_m(idxs) = antennas.raster.Coax_mm(idxs)*mm_to_m; % idxs that are NaN change to raster

idxs = find(isnan(output.Coax_m));
output.Coax_m(idxs) = antennas.install.Coax_mm(idxs)*mm_to_m; % idxs that are NaN change to install

idxs = find(isnan(output.Coax_m));
output.Coax_m(idxs) = antennas.ideal.Coax_mm(idxs)*mm_to_m; % idxs that are NaN change to ideal

% Now do bin parameters starting with height
if ~isnan(bin.best.height_ft)
    output.height_m = bin.best.height_ft*ft_to_m;
elseif ~isnan(bin.raster.height_ft)
    output.height_m = bin.raster.height_ft*ft_to_m;
elseif ~isnan(bin.install.height_ft)
    output.height_m = bin.install.height_ft*ft_to_m;
elseif isnan(bin.ideal.height_ft)
    output.height_m = bin.ideal.height_ft*ft_to_m;
else
    print('Could not determine bin height');
end

% Now do plenum
if ~isnan(bin.best.plenum_in)
    output.plenum_m = bin.best.plenum_in*in_to_m;
elseif ~isnan(bin.raster.plenum_in)
    output.plenum_m = bin.raster.plenum_in*in_to_m;
elseif ~isnan(bin.install.plenum_in)
    output.plenum_m = bin.install.plenum_in*in_to_m;
elseif isnan(bin.ideal.plenum_in)
    output.plenum_m = bin.ideal.plenum_in*in_to_m;
else
    print('Could not determine bin height');
end

% Now do diameter
if ~isnan(bin.best.diameter_ft)
    output.diameter_m = bin.best.diameter_ft*ft_to_m;
elseif ~isnan(bin.raster.diameter_ft)
    output.diameter_m = bin.raster.diameter_ft*ft_to_m;
elseif ~isnan(bin.install.diameter_ft)
    output.diameter_m = bin.install.diameter_ft*ft_to_m;
elseif isnan(bin.ideal.diameter_ft)
    output.diameter_m = bin.ideal.diameter_ft*ft_to_m;
else
    print('Could not determine diameter');
end

% Now do angle
if ~isnan(bin.best.angle_deg)
    output.angle_rad = bin.best.angle_deg*deg_to_rad;
elseif ~isnan(bin.raster.angle_deg)
    output.angle_rad = bin.raster.angle_deg*deg_to_rad;
elseif ~isnan(bin.install.angle_deg)
    output.angle_rad = bin.install.angle_deg*deg_to_rad;
elseif isnan(bin.ideal.angle_deg)
    output.angle_rad = bin.ideal.angle_deg*deg_to_rad;
else
    print('Could not angle');
end

% Now do cap
if ~isnan(bin.best.cap_diameter_ft)
    output.cap_diameter_m = bin.best.cap_diameter_ft*ft_to_m;
elseif ~isnan(bin.raster.cap_diameter_ft)
    output.cap_diameter_m = bin.raster.cap_diameter_ft*ft_to_m;
elseif ~isnan(bin.install.cap_diameter_ft)
    output.cap_diameter_m = bin.install.cap_diameter_ft*ft_to_m;
elseif isnan(bin.ideal.cap_diameter_ft)
    output.cap_diameter_m = bin.ideal.cap_diameter_ft*ft_to_m;
else
    print('Could not determine bin height');
end


%% Correct a problem with the way we compute angles

% We define the angle the antennas are at as the angle clockwise between
% the hatch and the antenna when in the bin (incrementing angle clockwise).
% We require a mathematical angle which increments counterclockwise
% We just use 2*pi - angle

output.Phi_rad = 2*pi - output.Phi_rad;


%% Generate output files


tx_fid = fopen(output_tx_file_name,'w');
rx_fid = fopen(output_rx_file_name,'w');
tx2rx_fid = fopen(output_tx_to_rx_file_name,'w');

fprintf(tx_fid,'%d\n',nTx);
fprintf(rx_fid,'%d\n',nTx);
fprintf(tx2rx_fid,'%d\n',nTx);


output.antennaX = (output.diameter_m/2 - 0.08).*cos(output.Phi_rad);
output.antennaY = (output.diameter_m/2 - 0.08).*sin(output.Phi_rad);
output.antennaZ = output.Z_m;
output.coordinates = [output.antennaX, output.antennaY output.antennaZ];

radius_to_antenna_vectors = [output.coordinates - [0, 0, 0]];
radius_to_antenna_vectors(:,3) = 0; % Since we are using this to to convert only to X-Y positions

Tx_J_vector = Tx_J_vector./norm(Tx_J_vector);
Rx_J_vector = Rx_J_vector./norm(Rx_J_vector);

disp('Moving antennas 8cm inward');
inward_shift = 0.08;
antenna_type = 30; % H-phi antennas
num_polarizations = 1;
rx_per_tx = nTx - 1;

for tx=1:24
    
    radius_to_antenna_vector = radius_to_antenna_vectors(tx,:)./norm(radius_to_antenna_vectors(tx,:));
    
    output.tx_polarizations(tx,[1:3]) = cross(Tx_J_vector, radius_to_antenna_vector);
    output.rx_polarizations(tx,[1:3]) = cross(Rx_J_vector, radius_to_antenna_vector);
    assert(abs(output.tx_polarizations(tx,3))<1e-6);
    
    
    fprintf(tx_fid,'%d\t%d\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\n',tx, antenna_type, ...
        output.coordinates(tx,1), output.coordinates(tx,2), output.coordinates(tx,3), ...
        output.tx_polarizations(tx,1),  output.tx_polarizations(tx,2),  output.tx_polarizations(tx,3));
    
    fprintf(rx_fid,'%d\t%1.8f\t%1.8f\t%1.8f\t%d\t%1.8f\t%1.8f\t%1.8f\n',tx, ...
        output.coordinates(tx,1), output.coordinates(tx,2), output.coordinates(tx,3), ...
        num_polarizations, ...
        output.rx_polarizations(tx,1),  output.rx_polarizations(tx,2),  output.rx_polarizations(tx,3));
    
    fprintf(tx2rx_fid,'%d\t%d\t',tx, rx_per_tx);
    for rx=1:24
        if tx ~= rx
            fprintf(tx2rx_fid,'%d\t',rx);
        end
    end
    fprintf(tx2rx_fid,'\n');
    
end

fclose(tx_fid);
fclose(rx_fid);
fclose(tx2rx_fid);

%% Read the files and plot the results to validate

Tx_pos = dlmread(output_tx_file_name,'\t',[1 2 24 4]);
Tx_pol = dlmread(output_tx_file_name,'\t',[1 5 24 7]);
nTx = size(Tx_pos,1);

Rx_pos = dlmread(output_rx_file_name,'\t',[1 1 24 3]); % this is different from Tx lines
Rx_pol = dlmread(output_rx_file_name,'\t',[1 5 24 7]);
nRx = size(Rx_pos,1);

figure(); hold on;

% Plot bin
eave_points = [0:0.01:2*pi].';
eave_points = [output.diameter_m/2*cos(eave_points), output.diameter_m/2*sin(eave_points), zeros(size(eave_points))];
plot3(eave_points(:,1), eave_points(:,2), eave_points(:,3), 'b-');

floor_points = eave_points - [0 0 output.height_m];
plot3(floor_points(:,1), floor_points(:,2), floor_points(:,3), 'b-');


% Plot antennas
tx_z_shift = -0.25;
rx_z_shift = -0.50;
pol_scale = 5;

for tx=1:nTx
    scatter3(Tx_pos(tx,1), Tx_pos(tx,2), Tx_pos(tx,3), 75, 'k', 'filled');
    text(Tx_pos(tx,1), Tx_pos(tx,2), Tx_pos(tx,3)+tx_z_shift, num2str(tx));
    h = quiver3(Tx_pos(tx,1),Tx_pos(tx,2), Tx_pos(tx,3), pol_scale*Tx_pol(tx,1),  pol_scale*Tx_pol(tx,2),  pol_scale*Tx_pol(tx,3), 'g');
end

for Rx=1:nRx
    scatter3(Rx_pos(Rx,1), Rx_pos(Rx,2), Rx_pos(Rx,3), 75, 'k', 'filled');
    text(Rx_pos(Rx,1), Rx_pos(Rx,2), Rx_pos(Rx,3)+rx_z_shift, num2str(Rx));
    g = quiver3(Rx_pos(Rx,1),Rx_pos(Rx,2), Rx_pos(Rx,3), pol_scale*Rx_pol(Rx,1),  pol_scale*Rx_pol(Rx,2),  pol_scale*Rx_pol(Rx,3),'r');
end
legend([h g],'Tx polarizations','Rx polarizations');
view(45,45);

