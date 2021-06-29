function [coax_lengths] = calculate_roof_coax(bin_diameter_ft, roof_angle, roof_fudge_factor, peak_to_tbox_ft, cable_run_diameter_ft, wire_cable_angles)
%{

calculates the length of coax to run from hanger bracket to top box

inputs:
diameter_ft - nominal bin diameter (not actual) [ft]
roof_angle - roof angle of bin [degrees]
roof_fudge_factor - multiplier for cable run from eave to peak
peak_to_tbox_ft - distance from roof peak to top box position [ft]
cable_run_diameter_ft - cable run diameter around roof cap [ft]

outputs:
coax_lengths - num_bins x num_antennas matrix, length of roof coax [mm]

%}

[num_bins,num_cables] = size(wire_cable_angles);
num_antennas = 2*num_cables;

odd_idx = 1:2:num_antennas-1;
even_idx = 2:2:num_antennas;


antenna_angles(:,odd_idx) = wire_cable_angles;
antenna_angles(:,even_idx) = wire_cable_angles;
antenna_angles(antenna_angles>pi) = 2*pi - antenna_angles(antenna_angles>pi); % restrict range to 0 to 180 degrees

roof_length = 0.3048*(bin_diameter_ft/2) ./ cosd(roof_angle); % length of roof

cable_up = roof_length.*roof_fudge_factor;
cable_down = 0.3048*peak_to_tbox_ft;

[m,~] = size(cable_run_diameter_ft);
if (m == 1) % make sure this is a column vector
    cable_run_diameter_ft = cable_run_diameter_ft*ones(num_bins,1); 
end

coax_lengths = 1000*(cable_up + repmat(cable_down,1,num_antennas) + (cable_run_diameter_ft*0.3048/2).*antenna_angles); % up + down + arc length based on cable angle and diameter

coax_lengths = 100*ceil(coax_lengths/100); % round to nearest 100 [mm]

end

