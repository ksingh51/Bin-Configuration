clear all;
% close all;
% clc;


diam_ft = 30; % nominal bin diameter in ft
roof_angle_rad = pi/6; % bin roof angle in radians
cable_run_diameter_ft = 5; % set diameter of cable run around roof cap
peak_to_tbox_ft = 12; % distance from bin peak to top box install location
bot_to_tbox_ft = 4.9; % distance from concrete base to top box install location
odd_idx = 1:2:23;
even_idx = 2:2:24;
cbl_angles_rad(odd_idx) = 0:pi/6:2*pi-pi/6;
cbl_angles_rad(even_idx) = cbl_angles_rad(odd_idx);
cbl_angles_rad(cbl_angles_rad>pi) = 2*pi - cbl_angles_rad(cbl_angles_rad>pi);

% TBoT

L1 = 0.3048*(diam_ft/2)/cos(roof_angle_rad);
L2 = 0.3048*(cable_run_diameter_ft/2)*cbl_angles_rad;
L3 = 0.3048*peak_to_tbox_ft;

Lt = L1 + L2 + L3;


% TBoB

L4 = 0.3048*(diam_ft/2)*cbl_angles_rad;
L5 = 0.3048*bot_to_tbox_ft;

Lb = L4 + L5;


% delta

Ld = Lb - Lt;

[Lt',Lb']