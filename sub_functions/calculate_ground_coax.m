function [coax_lengths] = calculate_ground_coax(tru_bin_diameter_ft, wire_cable_angles, corrugation, tb_angle, tb_height_mm)
%{

calculates the length of coax to run from bottom hanger bracket to top box

inputs:
bin_diameter_ft - nominal bin diameter (not actual) [ft]
wire_cable_angles - angular postion of each wire rope [rad]
corruagion - bin corrugation
tb_angle - angular position of top box mounting [rad]
tb_height_mm - install height of top box from bin floor (not plenum) [mm]

outputs:
coax_lengths - num_bins x num_antennas matrix, length of roof coax [mm]

%}

[~,num_cables] = size(wire_cable_angles);
num_antennas = 2*num_cables;

odd_idx = 1:2:num_antennas-1;
even_idx = 2:2:num_antennas;

% re-ref wire rope angles to top box angle
wire_cable_angles = wire_cable_angles - tb_angle;
wire_cable_angles(wire_cable_angles<0) = wire_cable_angles(wire_cable_angles<0) + 2*pi;

antenna_angles(:,odd_idx) = wire_cable_angles;
antenna_angles(:,even_idx) = wire_cable_angles;
antenna_angles(antenna_angles>pi) = 2*pi - antenna_angles(antenna_angles>pi); % restrict range to 0 to pi

tru_bin_radius_mm = tru_bin_diameter_ft*304.8/(2);

coax_arc = antenna_angles*tru_bin_radius_mm;

coax_lengths = coax_arc + tb_height_mm;

coax_lengths = 500*ceil(coax_lengths/500); % round to nearest 500 [mm]

end

