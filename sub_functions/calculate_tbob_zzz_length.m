function [zzz_lengths] = calculate_tbob_zzz_length(bin_eave_height_ft, bin_diameter_ft, roof_angle, roof_fudge_factor, tbox_fudge_length_ft)
%{

calculates the length of 6 conductor cable used in the system for TBoB

inputs:
bin_eave_height_ft - eave height of bin [ft]
diameter_ft - nominal bin diameter (not actual) [ft]
roof_angle - roof angle of bin [degrees]
roof_fudge_factor - multiplier for cable run from eave to peak
tbox_fudge_length_ft - length to add (+ve) or subtract (-ve) from cable
    length to account for install location (distance from base, distance 
    over from ladder, etc) [ft]

outputs:
zzz_lengths
.cbl_002049 - top to bottom box cable run [mm]
.cbl_002053 - temp cable extension [mm]

%}

roof_length = 0.3048*(bin_diameter_ft/2) ./ cosd(roof_angle); % length of roof
roof_cable = roof_length.*roof_fudge_factor;
wall_cable = bin_eave_height_ft*0.3048;
full_cable = roof_cable + wall_cable + 0.3048*tbox_fudge_length_ft;

% top to bottom box
zzz_lengths.cbl_002049 = 1000; % hard coded at 1000 mm for now
% temp cable extension
zzz_lengths.cbl_002053 = 250*ceil(1000*full_cable/250);

end

