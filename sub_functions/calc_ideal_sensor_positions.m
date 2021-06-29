function [positions] = calc_ideal_sensor_positions(eff_eave_height_ft, plenum_height_ft, floor_clearance_ft, eave_clearance_ft, num_antennas, antenna_length_ft)
%{

calculates the locations of the the wire rope ferrules on all cable strings

inputs:
eff_eave_height - effective eave height, lowest ring of bolt holes to highest ring of bolt holes [ft]
plenum_height_ft - distance from base of bin to plenum floor [ft]
diameter_ft - nominal bin diameter (not actual) [ft]
roof_angle - roof angle of bin [degrees]
floor_clearance_ft - distance from plenum to base of lowest antenna [ft]
eave_clearance_ft - distance from plenum to top of highest antenna [ft]
num_antennas - number of antennas, assumed 2 antennas per cable
antenna_length_ft - length of antenna body [ft]

outputs:
positions - num_bins x num_antennas matrix, distance from top of wire rope
loop to crimped ferrule [mm]

%}

% antenna and wire rope dimensions
antenna_length_m = antenna_length_ft*0.3048; % length of antenna body

height_m = (eff_eave_height_ft - plenum_height_ft)*0.3048; % distance from plenum floor to top ring of bolt holes

zmin = floor_clearance_ft*0.3048 + antenna_length_m/2; % z position of center of lowest antenna
zmax = height_m - eave_clearance_ft*0.3048 - antenna_length_m/2; % z position of center of highest antenna

step = (zmax - zmin)/(num_antennas - 1); % spacing between antennas

switch num_antennas
    case 24
        ant_index = 1:24;
        even_height_to_ant_idx = [1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24];
        ant_to_sensor_idx = [2, 1, 14, 13, 8, 7, 20, 19, 4, 3, 16, 15, 10, 9, 22, 21, 6, 5, 18, 17, 12, 11, 24, 23];
    case 32
        ant_index = 1:32;
        even_height_to_ant_idx = [1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23, 8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31, 16, 32];
        ant_to_sensor_idx = [2, 1, 18, 17, 10, 9, 26, 25, 4, 3, 20, 19, 12, 11, 28, 27, 6, 5, 22, 21, 14, 13, 30, 29, 8, 7, 24, 23, 16, 15, 32, 31];
end

even_heights = zmin + step*(ant_index-1);

antenna_heights(:,ant_index) = even_heights(:,even_height_to_ant_idx);

positions = round(1000*(height_m - antenna_heights(:,ant_to_sensor_idx))); % sensor center in mm, referenced to eave height


end

