clear all;
close all;
clc;

% behlen double checker

% AAC
eave_height_ft = 82.17;
plenum_height_ft = 12/12;
bin_diameter_ft = 106.94;
num_antennas = 24;
num_rings = 25;
antenna_length_ft = (2*356.35/1000)/0.3048; % center to ferrule is 356.35mm

bolt_to_brackethole_mm = 93.665;
brackethole_to_hanger_mm = 49.73;

sheet_arc_eff = 125.984/12; % (3.2m)
bolt_space = sheet_arc_eff/32; 
sheet_height_eff = 39.3701/12; % 1m


antpos_out = [1 	-11522    	0.006136	33700;     
2 	-22522    	0.006136	44700;     
3 	-5522     	0.527689	28200;     
4 	-16522    	0.527689	39200;     
5 	-8522     	1.055379	31700;     
6 	-19522    	1.055379	42700;     
7 	-2522     	1.576932	26200;     
8 	-13522    	1.576932	37200;     
9 	-10522    	2.098486	34700;     
10	-21522    	2.098486	45700;     
11	-4522     	2.626175	28700;     
12	-15522    	2.626175	39700;     
13	-7522     	3.147729	32200;    
14	-18522    	3.147729	43200;     
15	-1522     	3.663146	25700;     
16	-12522    	3.663146	36700;     
17	-9522     	4.196971	33700;     
18	-20522    	4.196971	44700;     
19	-3522     	4.718525	27200;     
20	-14522    	4.718525	38200;     
21	-6522     	5.233942	29700;     
22	-17522    	5.233942	40700;     
23	-1522     	5.767768	24200;     
24	-12522    	5.767768	35200];

true_to_eff_eave_delta = (eave_height_ft - sheet_height_eff*num_rings)/2;
sheet_centers_ft = [sheet_height_eff/2:sheet_height_eff:sheet_height_eff*num_rings - sheet_height_eff/2];
viable_antenna_centers_mm = round((sheet_centers_ft + true_to_eff_eave_delta)*304.8,0)';
viable_ferrule_positions_eave_mm = round((viable_antenna_centers_mm-(true_to_eff_eave_delta*304.8)-(antenna_length_ft*304.8/2)-bolt_to_brackethole_mm+brackethole_to_hanger_mm),0);

% AAJ
eave_height_ft = 92;
plenum_height_ft = 12/12;
bin_diameter_ft = 63.49;
num_antennas = 24;
num_rings = 28;
antenna_length_ft = (2*356.35/1000)/0.3048; % center to ferrule is 356.35mm

bolt_to_brackethole_mm = 93.665;
brackethole_to_hanger_mm = 49.73;

sheet_arc_eff = 125.984/12; % (3.2m)
bolt_space = sheet_arc_eff/32; 
sheet_height_eff = 39.3701/12; % 1m


antpos_out = [1 	-11522    	0.006136	33700;     
2 	-22522    	0.006136	44700;     
3 	-5522     	0.527689	28200;     
4 	-16522    	0.527689	39200;     
5 	-8522     	1.055379	31700;     
6 	-19522    	1.055379	42700;     
7 	-2522     	1.576932	26200;     
8 	-13522    	1.576932	37200;     
9 	-10522    	2.098486	34700;     
10	-21522    	2.098486	45700;     
11	-4522     	2.626175	28700;     
12	-15522    	2.626175	39700;     
13	-7522     	3.147729	32200;    
14	-18522    	3.147729	43200;     
15	-1522     	3.663146	25700;     
16	-12522    	3.663146	36700;     
17	-9522     	4.196971	33700;     
18	-20522    	4.196971	44700;     
19	-3522     	4.718525	27200;     
20	-14522    	4.718525	38200;     
21	-6522     	5.233942	29700;     
22	-17522    	5.233942	40700;     
23	-1522     	5.767768	24200;     
24	-12522    	5.767768	35200];

true_to_eff_eave_delta = (eave_height_ft - sheet_height_eff*num_rings)/2;
sheet_centers_ft = [sheet_height_eff/2:sheet_height_eff:sheet_height_eff*num_rings - sheet_height_eff/2];
viable_antenna_centers_mm = round((sheet_centers_ft + true_to_eff_eave_delta)*304.8,0)'
viable_ferrule_positions_eave_mm = round((viable_antenna_centers_mm-(true_to_eff_eave_delta*304.8)-(antenna_length_ft*304.8/2)-bolt_to_brackethole_mm+brackethole_to_hanger_mm),0)