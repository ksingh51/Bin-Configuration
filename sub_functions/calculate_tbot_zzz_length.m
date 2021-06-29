function [zzz_lengths] = calculate_tbot_zzz_length(bin_eave_height_ft, bin_diameter_ft, roof_angle, roof_fudge_factor, peak_to_tbox_ft)
%{

calculates the length of 6 conductor cable used in the system for TBoT

inputs:
bin_eave_height_ft - eave height of bin [ft]
diameter_ft - nominal bin diameter (not actual) [ft]
roof_angle - roof angle of bin [degrees]
roof_fudge_factor - multiplier for cable run from eave to peak
peak_to_tbox_ft - distance from roof peak to top box position [ft]

outputs:
zzz_lengths
.cbl_002049 - top to bottom box cable run [mm]
.cbl_002053 - temp cable extension [mm]
.cbl_002059 - peak sonar extension [mm]
.cbl_002060 - eave sonar extension [mm]

%}

roof_length = 0.3048*(bin_diameter_ft/2) ./ cosd(roof_angle); % length of roof

cable_up = roof_length.*roof_fudge_factor;
cable_down = 0.3048*peak_to_tbox_ft*ones(size(cable_up));

% top to bottom box
zzz_lengths.cbl_002049 = 250*ceil(1000*((0.3048*bin_eave_height_ft) + cable_up)/250);
% temp cable extension
zzz_lengths.cbl_002053 = 250*ceil(1000*cable_up/250);
% peak sonar extension
zzz_lengths.cbl_002059 = 250*ceil(1000*cable_down/250);
% eave sonar extension
zzz_lengths.cbl_002060 = 250*ceil(1000*cable_up/250);


end

