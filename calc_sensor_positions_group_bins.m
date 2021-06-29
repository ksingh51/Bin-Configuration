%--------------------------------------------------------------------------
%{
Sensor Position Calculator and Bin Grouper

151 Research Inc.
Kyle Nemez

This script:
1. Loads in bin information scraped from catalogues, formats it in a
struct, and fixes errors from scraping. Requires bindata.mat.
2. Basic plotting of bin sizes.
3. Groups bins according to number of rings, corrugation, and diameter
ranges.
4. Calculates vertical sensor positions. Ideal ones first, and then shifts
them so they fall in vertical ring centers.

The two outputs are structs: rbindata and gbindata

rbindata is raw bin data, contains all documented bin variants and their
physical data. specifically the subfields are:

rbindata
    .corrugation - categorical; Wide or Narrow or Behlen
    .manufacturer - categorical; GSI, WST, BRK, BLN
    .bottom - categorical; Flat or Hopper
    .num_rings - integer; number of vertical tiers/rings
    .eave_height_ft - double; the manufacturer specified eave height,
        assumed to be from floor to top of top tier [ft]
    .diameter_ft - double; nominal bin diameter listed by the manufacturer,
        note that the true diameter is slightly different [ft]
    .roof_angle - double; angle of the bin roof [degrees]
    .gvgsi_partno - string; unique GV GSI part number for the bin
    .antenna_length_ft - double; antenna center snip to ferrule clips
        multiplied by 2 [ft]
    .height_group - integer; the assigned height group from this script,
        can be used to map to gbindata
    .size_code - integer; the assigned size code from this script, can be
        used to map to gbindata
    .version - the version of this script that created the saved data

gbindata is grouped bin data. grouping is first done by height and then by
diameter. specifically the subfields are:

gbindata
    .size_code - integer; the assigned size code from this script, can be
        used to map to rbindata
    .min_gdiam_ft - double; minimum nominal diameter in the group [ft]
    .max_gdiam_ft - double; maximum nominal diameter in the group [ft]
    .height_group - integer; the assigned height group for the given size
        group, can be used to map to rbindata
    .num_rings - integer; number of vertical tiers/rings
    .eff_eave_height_ft - double; the effective height of the bin from the
        top ring of bolt holes to the bottom, found by multiplying the
        effective ring height (44" or 32") by the number of rings
    .corrugation - categorical; Wide or Narrow
    .antenna_length_ft - double; antenna center snip to ferrule clips
        multiplied by 2 [ft]
    .sensor_position_mm - double; row vector containing the 24 sensor
        positions for the size group. the distance is from the top ring of
        bolt holes (effective eave, not true eave) to the center of the
        antenna body. the antenna body center is equal to the center of the
        ring in which it resides [mm]
    .version - the version of this script that created the saved data

%}
%--------------------------------------------------------------------------
% Pleanum was always included and augur clearance was not included

tic;

clear all;
close all;
clc;

version = '1.2.5';

% script options
plot_general_sizes = 0;
plot_size_groups = 0;
compile_rawbindata = 1;
save_structs = 1;
simple_plot = 0;
use_manual_sensor_positions = 1;
evaluate_position_metrics = 0;

base_path = 'C:\Projects\MATLAB\GrainViz\151research-bin_configuration-dca7047f697a'; % get file path to script

%--------------------------------------------------------------------------
% ANY SIGNFICANT EDITS TO THIS SCRIPT SHOULD RESULT IN AN UPDATE TO THE
% VERSION NUMBER CONTAINED IN 'version' AND RE-SAVING OF THE STRUCTS
%--------------------------------------------------------------------------


%% load raw bin data and re-format

addpath(strcat(base_path,'\sub_functions'));

if (use_manual_sensor_positions)
    % load in pre-set sensor ring positions
    load(strcat(base_path,'\bin_structs\sensor_rings.mat')); %KD -  Updated by Kyle as of 6/1
end

if (~compile_rawbindata)
    % load pre-saved raw bin data
    load(strcat(base_path,'\bin_structs\rbindata.mat'));
    nn = length(rbindata.num_rings);
else
    load(strcat(base_path,'\raw_bin_data\bindata.mat'));
    
    % move tabular data into structs
    % manufacturer - Categorical, [GSI, WST, BRK, BLN]
    % corrugation - Categorical, [Wide, Narrow, Behlen]
    % bottom - Categorical, [Flat, Hopper]
    % num_rings - Integer
    % eave_height_ft - Double
    % diameter_ft - Double
    % roof_angle - Double
    
    % add missing bin (summit site) specs
    sm_corr = categorical({'Wide'});
    sm_bott = categorical({'Flat'});
    sm_nrng = 14;
    sm_ehtm = (51.4436)*0.3048;
    sm_dmft = 36;
    sm_rang = 30;
    
    % add missing bin (redwing site) specs
    rw_corr = categorical({'Narrow'});
    rw_bott = categorical({'Flat'});
    rw_nrng = 27;
    rw_ehtm = (72 + 2/12)*0.3048;
    rw_dmft = 156;
    rw_rang = 30;
    
    % add behlen bins at CFS
    cfs_corr = categorical({'Behlen','Behlen'}.');
    cfs_bott = categorical({'Flat';'Flat'});
    cfs_nrng = [25, 28].';
    cfs_ehtm = [(82+2/12)*0.3048, 92*0.3048].';
    cfs_dmft = [105, (62+4/12)].';
    cfs_rang = [29, 24.9].';
    
    bindata.corrugation = [bindatagsi.Corrugation;sm_corr;rw_corr;bindatawesteel.Corrugation;bindatabrock.Corrugation;bindatabrockstiff.Corrugation;cfs_corr(1);cfs_corr(2)];
    
    
    [n,~] = size(bindata.corrugation);
    [ngsi,~] = size(bindatagsi.Corrugation);
    ngsi = ngsi + 2; % +2 for redwing and summit
    [nwst,~] = size(bindatawesteel.Corrugation);
    [nbrk,~] = size([bindatabrock.Corrugation;bindatabrockstiff.Corrugation]);
    bindata.manufacturer = cell(n,1);
    bindata.manufacturer(1:ngsi,1) = {'GSI'};
    bindata.manufacturer((ngsi+1):(ngsi+nwst),1) = {'WST'};
    bindata.manufacturer((ngsi+nwst+1):(ngsi+nwst+nbrk),1) = {'BRK'};
    bindata.manufacturer((ngsi+nwst+nbrk+1):(n),1) = {'BLN'};
    bindata.manufacturer = categorical(bindata.manufacturer);
    
    % fix bad westeel data (scraping problems)
    wst_idx = find((bindatawesteel.Diameterft == 45) & (bindatawesteel.Corrugation == 'Wide'));
    correct_data = [13.0000   14.5800
        14.0000   15.7000
        15.0000   16.8200
        16.0000   17.9300
        17.0000   19.0500
        18.0000   20.1700
        19.0000   21.2900
        20.0000   22.4100
        21.0000   23.5300
        22.0000   24.6500
        23.0000   25.7600];
    bindatawesteel.EaveHeightm(wst_idx) = correct_data(:,2);
    
    wst_idx = find((bindatawesteel.Corrugation == 'Wide') & (bindatawesteel.Rings == 21) & (bindatawesteel.Diameterft == 102));
    bindatawesteel.EaveHeightm(wst_idx) = 23.53;
    
    wst_idx = find((bindatawesteel.Corrugation == 'Wide') & (bindatawesteel.Rings == 17) & (bindatawesteel.Diameterft == 108));
    bindatawesteel.EaveHeightm(wst_idx) = 19.06;
    
    wst_idx = find((bindatawesteel.Corrugation == 'Wide') & (bindatawesteel.Rings == 9) & (bindatawesteel.Diameterft == 84));
    bindatawesteel.EaveHeightm(wst_idx) = 10.11;
    
    wst_idx = find((bindatawesteel.Corrugation == 'Wide') & (bindatawesteel.Rings == 11) & (bindatawesteel.Diameterft == 105));
    bindatawesteel.EaveHeightm(wst_idx) = 12.35;
    
    wst_idx = find((bindatawesteel.Corrugation == 'Wide') & (bindatawesteel.Rings == 12) & (bindatawesteel.Diameterft == 105));
    bindatawesteel.EaveHeightm(wst_idx) = 13.46;
    
    bindata.bottom = [bindatagsi.Bottom;sm_bott;rw_bott;bindatawesteel.Bottom;bindatabrock.Bottom;bindatabrockstiff.Bottom;cfs_bott(1);cfs_bott(2)];
    bindata.num_rings = [bindatagsi.Rings;sm_nrng;rw_nrng;bindatawesteel.Rings;bindatabrock.Rings;bindatabrockstiff.Rings;cfs_nrng];
    bindata.eave_height_ft = [bindatagsi.EaveHeightm;sm_ehtm;rw_ehtm;bindatawesteel.EaveHeightm;bindatabrock.EaveHeightm;bindatabrockstiff.EaveHeightm;cfs_ehtm]/0.3048;
    bindata.diameter_ft = [bindatagsi.Diameterft;sm_dmft;rw_dmft;bindatawesteel.Diameterft;bindatabrock.Diameterft;bindatabrockstiff.Diameterft;cfs_dmft];
    bindata.roof_angle = [bindatagsi.RoofAngledegrees;sm_rang;rw_rang;bindatawesteel.RoofAngledegrees;bindatabrock.RoofAngledegrees;bindatabrockstiff.RoofAngledegrees;cfs_rang];
    
    % reduce data to flat bottom and minimum height
    min_height_ft = 20;
    good_bins_idx = find((bindata.bottom ~= 'Hopper') & ...
        ( (bindata.eave_height_ft >= min_height_ft) | ...
        (((bindata.eave_height_ft < 18.7) & (bindata.eave_height_ft > 18.4)) &  ((bindata.diameter_ft < 33) & (bindata.diameter_ft > 21))  ) ... % include extra bins for assumption illinois site
            ));
    nn = length(good_bins_idx);
    rbindata.corrugation = bindata.corrugation(good_bins_idx);
    rbindata.manufacturer = bindata.manufacturer(good_bins_idx);
    rbindata.bottom = bindata.bottom(good_bins_idx);
    rbindata.num_rings = bindata.num_rings(good_bins_idx);
    rbindata.eave_height_ft = bindata.eave_height_ft(good_bins_idx);
    rbindata.diameter_ft = bindata.diameter_ft(good_bins_idx);
    rbindata.roof_angle = bindata.roof_angle(good_bins_idx);
    
    % build GV GSI part number - KD - 6/1 - Needs updating but not sure if
    % I have the data reqd - Talk to kyle - Deerek will provide the data
    % for this
    for i = 1:length(rbindata.diameter_ft)
        if (rbindata.corrugation(i) == 'Wide')
            corr = 'W';
        elseif (rbindata.corrugation(i) == 'Narrow')
            corr = 'N';
        else
            corr = 'B';
        end
        
        if (isnan(rbindata.num_rings(i)))
            nr = '00';
        elseif (rbindata.num_rings(i) < 10)
            nr = strcat('0',num2str(rbindata.num_rings(i)));
        else
            nr = num2str(rbindata.num_rings(i));
        end
        
        if (rbindata.diameter_ft(i) < 100)
            dm = strcat('0',num2str(round(rbindata.diameter_ft(i),0)));
        else
            dm = num2str(round(rbindata.diameter_ft(i),0));
        end
        
        rbindata.gvgsi_partno{i,1} = strcat('GV',dm,corr,nr,char(rbindata.manufacturer(i)),'F','M');
    end
    
    clear bindatabrock bindatabrockstiff bindatagsi bindatawesteel bindata;
    
    % 
    
    % get index of nan rings
    nan_idx_binary = isnan(rbindata.num_rings);
    nan_idx = find(nan_idx_binary);
    
    % figure
    % hold on;
    % scatter(rbindata.diameter_ft(~nan_idx_binary),rbindata.eave_height_ft(~nan_idx_binary),'filled','b');
    % scatter(rbindata.diameter_ft(nan_idx_binary),rbindata.eave_height_ft(nan_idx_binary),'r');
    % xlabel('Diameter [ft]');
    % ylabel('Eave Height [ft]');
    % title('Bin Sizes');
    % legend('Ring Data','No Ring Data','Location','SouthEast');
    % xlim([0 160]);
    % ylim([0 110]);
    % grid on;
    
    % replace missing num rings in comparable bins (match for height within
    % 0.08', and same corrugation
    for i = 1:length(nan_idx)
        match_idx = find((abs(rbindata.eave_height_ft(~nan_idx_binary) - rbindata.eave_height_ft(nan_idx(i)))<0.1) & (rbindata.corrugation(~nan_idx_binary) == rbindata.corrugation(nan_idx(i))));
        if (~isempty(match_idx))
            rbindata.num_rings(nan_idx(i)) = rbindata.num_rings(match_idx(1));
        end
    end
    
    % find any remaining nan values
    nan_idx_binary = isnan(rbindata.num_rings);
    nan_idx = find(nan_idx_binary);
    
    % discovered that these are all Brock Wide Corr bins with no comparables
    % best guess for num rings is height/44" rounded down
    for i = 1:length(nan_idx)
        if (rbindata.corrugation(nan_idx(i)) == 'Wide')
            rings_guess = floor(rbindata.eave_height_ft(nan_idx(i))/(44/12));
            rbindata.num_rings(nan_idx(i)) = rings_guess;
        end
    end
    
    for i = 1:length(rbindata.num_rings)
        if (i == 1274)
            j = 2;
        end
        if (rbindata.num_rings(i) < 10)
            nr = strcat('0',num2str(rbindata.num_rings(i)));
        else
            nr = num2str(rbindata.num_rings(i));
        end
        rbindata.gvgsi_partno(i) = strrep(rbindata.gvgsi_partno(i),'00',nr);
    end
    
    % nan_idx_binary = isnan(rbindata.num_rings);
    % nan_idx = find(nan_idx_binary);
    %
    % figure
    % hold on;
    % scatter(rbindata.diameter_ft(~nan_idx_binary),rbindata.eave_height_ft(~nan_idx_binary),'filled','b');
    % scatter(rbindata.diameter_ft(nan_idx_binary),rbindata.eave_height_ft(nan_idx_binary),'r');
    % xlabel('Diameter [ft]');
    % ylabel('Eave Height [ft]');
    % title('Bin Sizes');
    % legend('Ring Data','No Ring Data','Location','SouthEast');
    % xlim([0 160]);
    % ylim([0 110]);
    % grid on;
end
% Change rbindata to delete bins with 12' dia - Update 6/20 (KD)
outdated_bins = rbindata.diameter_ft == 12 | rbindata.num_rings < 9; % bins with dia <15' and rings < 9
rbindata.corrugation(outdated_bins) = [];
rbindata.manufacturer(outdated_bins) = [];
rbindata.bottom(outdated_bins) = [];
rbindata.num_rings(outdated_bins) = [];
rbindata.eave_height_ft(outdated_bins) = [];
rbindata.roof_angle(outdated_bins) = [];
rbindata.gvgsi_partno(outdated_bins) = [];
rbindata.diameter_ft(outdated_bins) = [];
% FOR NOW change 132' to 135'
rbindata.diameter_ft(rbindata.diameter_ft == 132) = 135;
%% basic plotting

gsiw_idx = find((rbindata.corrugation == 'Wide') & (rbindata.manufacturer == 'GSI'));
gsin_idx = find((rbindata.corrugation == 'Narrow') & (rbindata.manufacturer == 'GSI'));
wstw_idx = find((rbindata.corrugation == 'Wide') & (rbindata.manufacturer == 'WST'));
wstn_idx = find((rbindata.corrugation == 'Narrow') & (rbindata.manufacturer == 'WST'));
brkw_idx = find((rbindata.corrugation == 'Wide') & (rbindata.manufacturer == 'BRK'));
brkn_idx = find((rbindata.corrugation == 'Narrow') & (rbindata.manufacturer == 'BRK'));
bln_idx = find((rbindata.corrugation == 'Behlen') & (rbindata.manufacturer == 'BLN'));

if plot_general_sizes
    figure
    hold on;
    scatter(rbindata.num_rings(gsiw_idx),rbindata.eave_height_ft(gsiw_idx),'filled','r');
    scatter(rbindata.num_rings(gsin_idx),rbindata.eave_height_ft(gsin_idx),'filled','r','d');
    scatter(rbindata.num_rings(wstw_idx),rbindata.eave_height_ft(wstw_idx),'filled','b');
    scatter(rbindata.num_rings(wstn_idx),rbindata.eave_height_ft(wstn_idx),'filled','b','d');
    scatter(rbindata.num_rings(brkw_idx),rbindata.eave_height_ft(brkw_idx),'filled','k');
    scatter(rbindata.num_rings(brkn_idx),rbindata.eave_height_ft(brkn_idx),'filled','k','d');
    scatter(rbindata.num_rings(bln_idx),rbindata.eave_height_ft(bln_idx),'filled','g');
    xlabel('Number of Rings');
    ylabel('Eave Height [ft]');
    title('Bin Heights');
    legend('GSI - Wide','GSI - Narrow','WST - Wide','WST - Narrow','BRK - Wide','BRK - Narrow','BLN - Behlen','Location','NorthWest');
    xlim([0 45]);
    ylim([0 110]);
    grid on;
    
    figure
    hold on;
    scatter(rbindata.diameter_ft(gsiw_idx),rbindata.eave_height_ft(gsiw_idx),'filled','r');
    scatter(rbindata.diameter_ft(gsin_idx),rbindata.eave_height_ft(gsin_idx),'filled','r','d');
    scatter(rbindata.diameter_ft(wstw_idx),rbindata.eave_height_ft(wstw_idx),'filled','b');
    scatter(rbindata.diameter_ft(wstn_idx),rbindata.eave_height_ft(wstn_idx),'filled','b','d');
    scatter(rbindata.diameter_ft(brkw_idx),rbindata.eave_height_ft(brkw_idx),'filled','k');
    scatter(rbindata.diameter_ft(brkn_idx),rbindata.eave_height_ft(brkn_idx),'filled','k','d');
    scatter(rbindata.diameter_ft(bln_idx),rbindata.eave_height_ft(bln_idx),'filled','g');
    xlabel('Diameter [ft]');
    ylabel('Eave Height [ft]');
    title('Bin Sizes');
    legend('GSI - Wide','GSI - Narrow','WST - Wide','WST - Narrow','BRK - Wide','BRK - Narrow','BLN - Behlen','Location','SouthEast');
    xlim([0 160]);
    ylim([0 110]);
    grid on;
end

%% define parameters and constraints

plenum_height_ft = 18/12; %(UPDATE 6/1 KD - changed plenum ht to 18). worst case plenum flor height is 17.25", this is referenced to lowest ring of bolt holes on bottom sheet
num_antennas = 24; % number of antennas is 24
num_cables = num_antennas/2; % number of wire rope cables assemblies is 12

w_sheet_height_ft = 44/12; % wide corr sheet height is 44"
n_sheet_height_ft = 32/12; % narrow corr sheet height is 32"
b_sheet_height_ft = 39.3701/12; % behlen corr sheet height is 39-1/4" (1m)

eave_clearance_ft = 32/12; % clearance from eave to highest antenna top, now one narrow ring, used to be 0.4826/0.3048
floor_clearance_ft = 0.762/0.3048; % clearance from plenum floor to lowest  bottom

s_ant_length = (2*211.85/1000)/0.3048; % center to ferrule is 211.85mm
m_ant_length = (2*356.35/1000)/0.3048; % center to ferrule is 356.35mm
l_ant_length = (2*495.35/1000)/0.3048; % center to ferrule is 495.35mm

% 6/17 - Effective eave height calculated for all bins instead of groups
rbindata.eff_eave_height_ft = zeros(length(rbindata.corrugation),1);
rbindata.eff_eave_height_ft(rbindata.corrugation == 'Wide',1) = rbindata.num_rings(rbindata.corrugation == 'Wide').*w_sheet_height_ft;
rbindata.eff_eave_height_ft(rbindata.corrugation == 'Narrow',1) = rbindata.num_rings(rbindata.corrugation == 'Narrow').*n_sheet_height_ft;
rbindata.eff_eave_height_ft(rbindata.corrugation == 'Behlen',1) = rbindata.num_rings(rbindata.corrugation == 'Behlen').*b_sheet_height_ft;

%% antenna length ft based on bin size and dia (Updated 6/18 - KD:  details in Sensor Size matrix sheet)
% bin diameter <= 50, small antenna
% bin diameter >50 and <=90, medium antenna
% bin diameter >90, large antenna
% except for behlen bins, which get either small or medium, never large

sensor_size_ma
trix = {'',9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40;
    15,'Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','','','','','','','','','','','','','','','','','','';
    18,'Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','','','','','','','','','','','','','','','','','','';
    21,'Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','','','','','','','','','','','','','','','','';
    24,'Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','','','','','','','','','','','','';
    27,'Small','Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','','','','','','','','','','';
    30,'Small','Small','Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','','','','','','','','';
    33,'Small','Small','Small','Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','','','','','','','','';
    36,'Small','Small','Small','Small','Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','','','','','','','','';
    42,'Small','Small','Small','Small','Small','Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','','','','','','';
    48,'Small','Small','Small','Small','Small','Small','Small','Small','Small','Small','Small','Small','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','','','','','','';
    54,'Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','','','';
    60,'Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large';
    72,'Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large';
    75,'Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large';
    78,'Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large';
    90,'Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Medium','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large';
    105,'Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large';
    135,'Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','','','','','','';
    156,'Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','Large','','','','','','','','','','','','','';};

rbindata.antenna_length_ft = NaN*zeros(length(rbindata.corrugation),1);

% adding values for bottom sensor override. This is a test using data for
% 3/16 array with plenum  (in inches)
bottom_clearance_sm = 28.36;
bottom_clearance_med = 38.86;
bottom_clearance_lg = 33.36;


for i = 1:length(rbindata.corrugation)
   if rbindata.corrugation(i,1) == 'Narrow'
       if strcmp(sensor_size_matrix{find([sensor_size_matrix{:,1}] == rbindata.diameter_ft(i,1)),...
               find([sensor_size_matrix{1,:}] == rbindata.num_rings(i,1))},'Small')
           rbindata.antenna_length_ft(i,1) = s_ant_length;
       elseif strcmp(sensor_size_matrix{find([sensor_size_matrix{:,1}] == rbindata.diameter_ft(i,1)),...
               find([sensor_size_matrix{1,:}] == rbindata.num_rings(i,1))},'Medium')
           rbindata.antenna_length_ft(i,1) = m_ant_length;
       else
           rbindata.antenna_length_ft(i,1) = l_ant_length;
       end
   elseif rbindata.corrugation(i,1) == 'Wide'
       rbindata.antenna_length_ft(i,1) = m_ant_length;

   end
       
           
end

% % assign antenna size to bin parameters
% rbindata.antenna_length_ft = NaN*zeros(nn,1);
% rbindata.antenna_length_ft(s_ant_idx) = s_ant_length;
% rbindata.antenna_length_ft(m_ant_idx) = m_ant_length;
% rbindata.antenna_length_ft(l_ant_idx) = l_ant_length;

%% height grouping

% find all unique numbers of rings for wide corr bins
kw = unique(rbindata.num_rings(rbindata.corrugation == 'Wide'));
kw = kw(~isnan(kw));
num_kw = length(kw);

% find all unique numbers of rings for narrow corr bins
kn = unique(rbindata.num_rings(rbindata.corrugation == 'Narrow'));
kn = kn(~isnan(kn));
num_kn = length(kn);

% find all unique numbers of rings for behlen corr bins
kb = unique(rbindata.num_rings(rbindata.corrugation == 'Behlen'));
kb = kb(~isnan(kb));
num_kb = length(kb);

num_hgroups = num_kw + num_kn + num_kb;

% start building new struct containing height group info
hbindata.num_rings = zeros(num_hgroups,1);
hbindata.height_group = zeros(num_hgroups,1);

hbindata.num_rings = [kw;kn;kb];
hbindata.height_group = [1:1:num_hgroups]'; % assign height groups, wide then narrow
hbindata.corrugation = cell(num_hgroups,1);
hbindata.corrugation(1:num_kw,1) = {'Wide'};
hbindata.corrugation(num_kw+1:num_kw+num_kn,1) = {'Narrow'};
hbindata.corrugation(num_kw+num_kn+1:num_kw+num_kn+num_kb,1) = {'Behlen'};
hbindata.corrugation = categorical(hbindata.corrugation);

rbindata.height_group = zeros(nn,1);

% establish height group mapping between rbindata & hbindata
for i = 1:num_kw
    ii = (rbindata.num_rings == kw(i)) & (rbindata.corrugation == 'Wide');
    rbindata.height_group(ii,1) = i;
end
for j = num_kw+1:num_kw+num_kn
    jj = (rbindata.num_rings == kn(j-num_kw)) & (rbindata.corrugation == 'Narrow');
    rbindata.height_group(jj,1) = j;
end
for k = num_kw+num_kn+1:num_kw+num_kn+num_kb
    kk = (rbindata.num_rings == kb(k-(num_kw+num_kn))) & (rbindata.corrugation == 'Behlen');
    rbindata.height_group(kk,1) = k;
end


% determine effective eave height (top bolts holes to bottom bolt holes)
hbindata.eff_eave_height_ft(1:num_kw) = hbindata.num_rings(1:num_kw).*w_sheet_height_ft;
hbindata.eff_eave_height_ft(1+num_kw:num_kw+num_kn) = hbindata.num_rings(1+num_kw:num_kw+num_kn).*n_sheet_height_ft;
hbindata.eff_eave_height_ft(1+num_kw+num_kn:num_kw+num_kn+num_kb) = hbindata.num_rings(1+num_kw+num_kn:num_kw+num_kn+num_kb).*b_sheet_height_ft;

%% Commented out - bin size grouping
% %% bin size grouping - THIS WILL MOST PROB GO AWAY
% 
% diameter_group_range = 21; % size groups are 21 ft
% bin_size_group = zeros(nn,1);
% last_max = 0;
% for i = 1:num_hgroups
%     idx = find(rbindata.height_group == i);
%     % get diameter group edges for the given height group
%     min_d = min(rbindata.diameter_ft(idx));
%     max_d = max(rbindata.diameter_ft(idx));
%     
%     % force edge to comply with antenna size bounds
%     % prevents from two bins in the same size group from getting two
%     % different antenna sizes
%     if (min_d < sm_break)
%         if (max_d <= sm_break)
%             diameter_edges = [min_d:diameter_group_range:max_d,max_d];
%         elseif ((max_d > sm_break) && (max_d <= ml_break))
%             diameter_edges = [min_d:diameter_group_range:sm_break,sm_break:diameter_group_range:max_d,max_d];
%         else % max d is strictly greater than the m/l break size
%             diameter_edges = [min_d:diameter_group_range:sm_break,sm_break:diameter_group_range:ml_break,ml_break:diameter_group_range:max_d,max_d];
%         end
%     elseif (min_d >= sm_break)
%         if (max_d <= sm_break)
%             diameter_edges = [min_d:diameter_group_range:max_d,max_d];
%         elseif ((max_d > sm_break) && (max_d <= ml_break))
%             diameter_edges = [min_d:diameter_group_range:max_d,max_d];
%         else % max d is strictly greater than the m/l break size
%             diameter_edges = [min_d:diameter_group_range:ml_break,ml_break:diameter_group_range:max_d,max_d];
%         end
%     elseif (min_d >= ml_break)
%         diameter_edges = [min_d:diameter_group_range:max_d,max_d];
%     end
%     
%     % sort data into bins based on diameter_edges
%     bin_size_group(idx) = discretize(rbindata.diameter_ft(idx),diameter_edges,'IncludedEdge','right'); % get diameter size group numbers
%     
%     % check if there are some diameter sizes that don't exist
%     if (length(unique(bin_size_group(idx))) ~= (length(diameter_edges)-1))
%         for j = 1:length(diameter_edges)-1
%             if (isempty(find(bin_size_group(idx)==j))) % find which one is missing
%                 idx_to_shift = find(bin_size_group(idx)>j);
%                 bin_size_group(idx(idx_to_shift)) = bin_size_group(idx(idx_to_shift)) - 1; % decrement diameter size
%                 j = j - 1; % decrement counter just in case this index is still missing
%             end
%         end
%     end
%     
%     % check for cases where there is a gap in sizes such that there are two
%     % neighbouring diameter groups where the min of the lower one and the
%     % max of the upper one are less than diameter_group_range apart. in
%     % that case, combine those groups, and decrement all higher groups by 1
%     if (max(bin_size_group(idx))>=2)
%         for j = 1:max(bin_size_group(idx))-1
%             
%             % find indexes of bins in the lower diam group
%             sub_idx = idx(bin_size_group(idx) == j);
%             lo_diams = rbindata.diameter_ft(sub_idx);
%             % get the smallest diameter bin in that group
%             loest_diam = min(lo_diams);
%             
%             % find indexes of bins in the upper diam group
%             sub_idx = idx(bin_size_group(idx) == j+1);
%             hi_diams = rbindata.diameter_ft(sub_idx);
%             % get the largest diameter bin in that group
%             hiest_diam = max(hi_diams);
%             
%             % if the high low difference is less than group range
%             if ((hiest_diam-loest_diam)<diameter_group_range)
%                 % make sure the collapsing groups don't collapse over sm/ml
%                 % antenna size breaks
%                 if ~(((hiest_diam > sm_break) && (loest_diam <= sm_break)) || ((hiest_diam > ml_break) && (loest_diam <= ml_break)))
%                     % find all bins in higher group numbers
%                     sub_idx = idx(bin_size_group(idx) >= j+1);-
%                     % decrement
%                     bin_size_group(sub_idx) = bin_size_group(sub_idx) - 1;
%                 end
%             end
%         end
%     end
%     
%     bin_size_group(idx) = bin_size_group(idx) + last_max; % add last maximum bin size value to ensure continuous distribution of bin sizes
%     last_max = max(bin_size_group(idx)); % update last max
%     
% end
% 
% num_bin_sizes = max(bin_size_group);
% rbindata.size_code = bin_size_group;
% 
% % build struct for size groups
% gbindata.size_code = zeros(num_bin_sizes,1);
% gbindata.size_code = [1:num_bin_sizes]';
% gbindata.min_gdiam_ft = zeros(num_bin_sizes,1);
% gbindata.max_gdiam_ft = zeros(num_bin_sizes,1);
% gbindata.height_group = zeros(num_bin_sizes,1);
% 
% for i = 1:num_bin_sizes
%     gbindata.height_group(i) = mean(rbindata.height_group(rbindata.size_code == i)); % these should all be the same
%     gbindata.min_gdiam_ft(i) = min(rbindata.diameter_ft(rbindata.size_code == i));
%     gbindata.max_gdiam_ft(i) = max(rbindata.diameter_ft(rbindata.size_code == i));
% end
% 
% gbindata.num_rings = zeros(num_bin_sizes,1);
% gbindata.eff_eave_height_ft = zeros(num_bin_sizes,1);
% gbindata.corrugation = cell(num_bin_sizes,1);
% 
% for i = 1:num_hgroups
%     hh = (gbindata.height_group==i);
%     gbindata.num_rings(hh) = hbindata.num_rings(i);
%     gbindata.eff_eave_height_ft(hh) = hbindata.eff_eave_height_ft(i);
%     gbindata.corrugation(hh) = cellstr(hbindata.corrugation(i));
% end
% 
% gbindata.corrugation = categorical(gbindata.corrugation);
% 
% % get size codes of behlen bins
% temp = zeros(num_bin_sizes,1);
% temp2 = rbindata.size_code(rbindata.manufacturer=='BLN');
% temp(temp2,1) = 1;
% 
% s_ant_idx = gbindata.max_gdiam_ft<=sm_break;
% m_ant_idx = (((gbindata.min_gdiam_ft>sm_break) & (gbindata.max_gdiam_ft<=ml_break)) | ((gbindata.min_gdiam_ft>sm_break) & temp));
% l_ant_idx = (gbindata.min_gdiam_ft>ml_break) & ~temp;
% 
% % assign antenna size to bin parameters
% gbindata.antenna_length_ft = zeros(num_bin_sizes,1);
% gbindata.antenna_length_ft(s_ant_idx) = s_ant_length;
% gbindata.antenna_length_ft(m_ant_idx) = m_ant_length;
% gbindata.antenna_length_ft(l_ant_idx) = l_ant_length;
% if (gbindata.antenna_length_ft == 0)
%     print('something is bust');
% end
% 
% if (plot_size_groups)
%     % PLOTTING
%     for i = 1:num_bin_sizes
%         wcbin_size_idx{i} = find((rbindata.size_code == i) & (rbindata.corrugation == 'Wide'));
%         ncbin_size_idx{i} = find((rbindata.size_code == i) & (rbindata.corrugation == 'Narrow'));
%     end
%     
%     color_loop_value = 10;
%     plot_cols = lines(num_bin_sizes);
%     cidx = 1;
%     cseed = 1;
%     figure('Position',[100 100 1500 1000]);
%     subplot(1,2,1);
%     hold on
%     for i = 1:num_bin_sizes
%         scatter(rbindata.diameter_ft(wcbin_size_idx{i}),rbindata.eave_height_ft(wcbin_size_idx{i}),50,'filled','MarkerFaceColor',plot_cols(cidx,:),'MarkerEdgeColor','k')
%         cidx = cidx + color_loop_value;
%         if (cidx >= length(plot_cols))
%             cidx = cseed;
%             cseed = cseed + 1;
%         end
%     end
%     xlabel('Diameter [ft]');
%     ylabel('Eave Height [ft]');
%     title('Wide Corr Bin Sizes');
%     xlim([0 160]);
%     ylim([0 110]);
%     grid on;
%     
%     subplot(1,2,2);
%     hold on
%     for i = 1:num_bin_sizes
%         scatter(rbindata.diameter_ft(ncbin_size_idx{i}),rbindata.eave_height_ft(ncbin_size_idx{i}),50,'filled','MarkerFaceColor',plot_cols(cidx,:),'MarkerEdgeColor','k')
%         cidx = cidx + color_loop_value;
%         if (cidx >= length(plot_cols))
%             cidx = cseed;
%             cseed = cseed + 1;
%         end
%     end
%     xlabel('Diameter [ft]');
%     ylabel('Eave Height [ft]');
%     title('Narrow Corr Bin Sizes');
%     xlim([0 160]);
%     ylim([0 110]);
%     grid on;
% end
%% Commented out  - sensor positions grouped
% %% calculate sensor positions - grouped % UPDATE 6/15 - commented away
% 
% % start with ideal positions, referenced from top ring of bolt holes
% %ideal_positions_mm = calc_ideal_sensor_positions(gbindata.eff_eave_height_ft,plenum_height_ft,floor_clearance_ft,eave_clearance_ft,num_antennas, gbindata.antenna_length_ft);
% num_bins = length(rbindata.antenna_length_ft);
% % define parameters
% stiff_per_sheet = 3;
% alt_stagger = 1;
% plot_2d = 0;
% plot_3d = 0;
% 
% plot_cols = parula(10);
% % UPDAT KD 6/14. No more grouping so sensor position array changes
% %sensor_position_ft = zeros(num_bin_sizes,num_antennas); 
% sensor_position_ft = zeros(num_bins, num_antennas);
% if (use_manual_sensor_positions)
%     for i = 1:num_bin_sizes
%         % get ring centers, referenced from bottom ring of bolts
%         if (gbindata.corrugation(i) == 'Wide')
%             ring_height = w_sheet_height_ft;
%             idx = sensor_rings.wide.ring_id(:,(sensor_rings.wide.nrings==gbindata.num_rings(i)));
%         elseif (gbindata.corrugation(i) == 'Narrow')
%             ring_height = n_sheet_height_ft;
%             idx = sensor_rings.narrow.ring_id(:,(sensor_rings.narrow.nrings==gbindata.num_rings(i)));
%         else
%             ring_height = b_sheet_height_ft;
%             idx = sensor_rings.behlen.ring_id(:,(sensor_rings.behlen.nrings==gbindata.num_rings(i)));
%         end
%         ring_centers_ft = gbindata.num_rings(i)*ring_height - ring_height/2:-1*ring_height:ring_height/2;%ring_height/2:ring_height:gbindata.num_rings(i)*ring_height;
%         
%         if (gbindata.num_rings(i) == 5)
%             % get half ring locations
%             ring_centers_mod_ft = [ring_centers_ft(1), mean([ring_centers_ft(1),ring_centers_ft(2)]), ring_centers_ft(2), mean([ring_centers_ft(2),ring_centers_ft(3)]), ring_centers_ft(3), mean([ring_centers_ft(3),ring_centers_ft(4)]), ring_centers_ft(4), mean([ring_centers_ft(4),ring_centers_ft(5)]), ring_centers_ft(5)];
%             jdx = [4,7,3,6,4,7,3,6,4,7,2,5,3,6,2,5,4,7,2,5,3,6,2,5];
%             sensor_position_ft(i,:) = ring_centers_mod_ft(jdx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
%             sensor_position_ft(i,[15,23]) = sensor_position_ft(i,[15,23]) + (16/12);
%         else
%             sensor_position_ft(i,:) = ring_centers_ft(idx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
%         end
%        
% %         sensor_position_ft(i,:) = ring_centers_ft(idx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
%         gbindata.sensor_position_mm(i,:) = (gbindata.eff_eave_height_ft(i)-sensor_position_ft(i,:))*304.8; % reref to top bolt hole, mm
%     end
%     
% else
%     for i = 1:num_bin_sizes
%         
%         % get ring centers, referenced from bottom ring of bolts
%         if (gbindata.corrugation(i) == 'Wide')
%             ring_height = w_sheet_height_ft;
%         elseif (gbindata.corrugation(i) == 'Narrow')
%             ring_height = n_sheet_height_ft;
%         else
%             ring_height = b_sheet_height_ft;
%         end
%         ring_centers_ft = ring_height/2:ring_height:gbindata.num_rings(i)*ring_height;
%         
%         % only accept ring centers that are greater than the plenum floor +
%         % clearance + half an antenna body
%         ring_centers_ft = ring_centers_ft(ring_centers_ft>=(plenum_height_ft + floor_clearance_ft + gbindata.antenna_length_ft(i)/2));
%         ring_centers_ft = ring_centers_ft(ring_centers_ft<=(gbindata.eff_eave_height_ft(i) - eave_clearance_ft - gbindata.antenna_length_ft(i)/2));
%         
%         % shift to nearest sheet center
%         reref_antenna_positions_ft = gbindata.eff_eave_height_ft(i) - ideal_positions_mm(i,:)/304.8; % rerefence the ideal antenna positions to ft from bottom bolt
%         
%         if (plot_2d)
%             [h_2d,h_3d,tru_diam_ft,~,theta_bolt,~] = plot_bin_panels(gbindata.num_rings(i),gbindata.corrugation(i),gbindata.max_gdiam_ft(i),plenum_height_ft,3,1,1,0);
%             
%             th_start = theta_bolt(1,1);
%             th_ideal_step = 2*pi/(num_antennas/2); % ideal distribution is even 30deg spacing (pi/6), referenced to first bolt hole for 24 antennas, 22.5 deg for 32 antennas
%             th_cables = th_start:th_ideal_step:th_start + 2*pi - th_ideal_step;
%             r_cables = (tru_diam_ft/2)*ones(size(th_cables)) - 4/12; % shift in by 4"
%             [x_cables, y_cables] = pol2cart(th_cables, r_cables);
%             xy_vals = th_cables.*(tru_diam_ft/2);
%             
%             set(0,'CurrentFigure',h_2d);
%             r = tru_diam_ft/2;
%             conv = 180/(pi*r);
%             
% %             for j = 1:num_cables
% %                 plot([xy_vals(j)*conv,xy_vals(j)*conv],[gbindata.eff_eave_height_ft(i), gbindata.eff_eave_height_ft(i) - gbindata.eff_eave_to_wirerope_bottom_mm(i)/304.8],'LineWidth',2,'Color',plot_cols(1,:));
% %             end
%             scatter(xy_vals*conv,reref_antenna_positions_ft(1:2:num_antennas-1),120,'Filled','s','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
%             scatter(xy_vals*conv,reref_antenna_positions_ft(2:2:num_antennas),120,'Filled','d','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
%             
%             %
%             %         t1 = repmat(theta_bolt*tru_diam_ft/2,length(ring_centers_ft),1);
%             %         m = length(theta_bolt);
%             %         t2 = repmat(ring_centers_ft,1,m);
%             %         scatter(t1(:)*conv,t2(:),20,'d','r','Filled')
%             %         title(rbindata.gvgsi_partno{i});
%             %
%             
%         end
%         
%         ip_mtx = repmat(reref_antenna_positions_ft',1,length(ring_centers_ft)); % make into matrix
%         rc_mtx = repmat(ring_centers_ft,length(reref_antenna_positions_ft'),1);
%         delta = abs(ip_mtx - rc_mtx);
%         [mins,idx] = min(delta,[],2); % find the closest ring center
%         
%         idx(idx==0) = idx(idx==0)+1;
%         idx(idx>length(ring_centers_ft)) = idx(idx>length(ring_centers_ft))-1;
%         
%         sensor_position_ft(i,:) = ring_centers_ft(idx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
%         gbindata.sensor_position_mm(i,:) = (gbindata.eff_eave_height_ft(i)-sensor_position_ft(i,:))*304.8; % reref to top bolt hole, mm
%     end
% end
% 
% odd_idx = [1:2:23];
% even_idx = [2:2:24];
% 
% % plot all distributions
% for j = 1:num_hgroups
%     gdx = find(gbindata.height_group==j);
%     i = gdx(1);
%     
%     if (simple_plot)
%         
%         if (gbindata.corrugation(i) == 'Wide')
%             psh_ft = w_sheet_height_ft;
%             ttle = 'Wide';
%         elseif (gbindata.corrugation(i) == 'Narrow')
%             psh_ft = n_sheet_height_ft;
%             ttle = 'Narrow';
%         else
%             psh_ft = b_sheet_height_ft;
%             ttle = 'Behlen';
%         end
%         plot_cols = parula(10);
%         
%         figure;
%         hold on;
%         line([0 13],[0, 0],'Color','blue','LineWidth',1.5);
%         line([0 13],[plenum_height_ft, plenum_height_ft],'Color','blue','LineWidth',1.5);
%         line([0 13],[gbindata.eff_eave_height_ft(i),gbindata.eff_eave_height_ft(i)],'Color','blue','LineWidth',1.5);
%         for k = 1:gbindata.num_rings(i)
%             line([0 13],[psh_ft*k,psh_ft*k],'Color','black');
%         end
%         
%         scatter([1:1:12],sensor_position_ft(i,odd_idx),'Filled','s','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
%         scatter([1:1:12],sensor_position_ft(i,even_idx),'Filled','d','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
%         
%         % grid on;
%         xlim([0 13]);
%         xticks([1:1:12])
%         yticks([psh_ft/2:psh_ft:(gbindata.num_rings(i)*psh_ft - psh_ft/2)])
%         ylim([gbindata.eff_eave_height_ft(i)*-0.1 gbindata.eff_eave_height_ft(i)*1.1])
%         xlabel('Cable Number');
%         ylabel('Ring ID [ft]');
% %         title(strcat(ttle," Corrugation, ", num2str(gbindata.num_rings(i))," Rings"));
%         
%         yt = gbindata.num_rings(i):-1:1;
%         yticklabels(compose('%d',yt));
%         
%     end
%     
% end

%% calculate sensor positions - ungrouped bins

% start with ideal positions, referenced from top ring of bolt holes
%ideal_positions_mm = calc_ideal_sensor_positions(gbindata.eff_eave_height_ft,plenum_height_ft,floor_clearance_ft,eave_clearance_ft,num_antennas, gbindata.antenna_length_ft);
num_bins = length(rbindata.antenna_length_ft);
% define parameters
stiff_per_sheet = 3;
alt_stagger = 1;
plot_2d = 0;
plot_3d = 0;

plot_cols = parula(10);
% UPDAT KD 6/14. No more grouping so sensor position matrix changes
%sensor_position_ft = zeros(num_bin_sizes,num_antennas); 
sensor_position_ft = zeros(num_bins, num_antennas);
if (use_manual_sensor_positions)
    for i = 1:num_bins
        % get ring centers, referenced from bottom ring of bolts
        if (rbindata.corrugation(i) == 'Wide')
            ring_height = w_sheet_height_ft;
            idx = sensor_rings.wide.ring_id(:,(sensor_rings.wide.nrings==rbindata.num_rings(i)));
        elseif (rbindata.corrugation(i) == 'Narrow')
            ring_height = n_sheet_height_ft;
            idx = sensor_rings.narrow.ring_id(:,(sensor_rings.narrow.nrings==rbindata.num_rings(i)));
        else
            ring_height = b_sheet_height_ft;
            idx = sensor_rings.behlen.ring_id(:,(sensor_rings.behlen.nrings==rbindata.num_rings(i)));
        end
        ring_centers_ft = rbindata.num_rings(i)*ring_height - ring_height/2:-1*ring_height:ring_height/2;%ring_height/2:ring_height:gbindata.num_rings(i)*ring_height;
        
        if (rbindata.num_rings(i) == 5) % dont change anything or just remove it
            % get half ring locations
            ring_centers_mod_ft = [ring_centers_ft(1), mean([ring_centers_ft(1),ring_centers_ft(2)]), ring_centers_ft(2), mean([ring_centers_ft(2),ring_centers_ft(3)]), ring_centers_ft(3), mean([ring_centers_ft(3),ring_centers_ft(4)]), ring_centers_ft(4), mean([ring_centers_ft(4),ring_centers_ft(5)]), ring_centers_ft(5)];
            jdx = [4,7,3,6,4,7,3,6,4,7,2,5,3,6,2,5,4,7,2,5,3,6,2,5];
            sensor_position_ft(i,:) = ring_centers_mod_ft(jdx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
            sensor_position_ft(i,[15,23]) = sensor_position_ft(i,[15,23]) + (16/12);
        else
            sensor_position_ft(i,:) = ring_centers_ft(idx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
        end
       
%         sensor_position_ft(i,:) = ring_centers_ft(idx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
         gbindata.sensor_position_mm(i,:) = (rbindata.eff_eave_height_ft(i)-sensor_position_ft(i,:))*304.8; % reref to top bolt hole, mm
    end  
else
    for i = 1:num_bins
 
        % get ring centers, referenced from bottom ring of bolts
        if (rbindata.corrugation(i) == 'Wide')
            ring_height = w_sheet_height_ft;
        elseif (rbindata.corrugation(i) == 'Narrow')
            ring_height = n_sheet_height_ft;
        else
            ring_height = b_sheet_height_ft;
        end
        ring_centers_ft = ring_height/2:ring_height:rbindata.num_rings(i)*ring_height;
        
        % only accept ring centers that are greater than the plenum floor +
        % clearance + half an antenna body
        ring_centers_ft = ring_centers_ft(ring_centers_ft>=(plenum_height_ft + floor_clearance_ft + rbindata.antenna_length_ft(i)/2));
        ring_centers_ft = ring_centers_ft(ring_centers_ft<=(rbindata.eff_eave_height_ft(i) - eave_clearance_ft - rbindata.antenna_length_ft(i)/2));
        
        % shift to nearest sheet center
        reref_antenna_positions_ft = rbindata.eff_eave_height_ft(i) - ideal_positions_mm(i,:)/304.8; % rerefence the ideal antenna positions to ft from bottom bolt
        
        if (plot_2d)
            [h_2d,h_3d,tru_diam_ft,~,theta_bolt,~] = plot_bin_panels(rbindata.num_rings(i),rbindata.corrugation(i),gbindata.max_gdiam_ft(i),plenum_height_ft,3,1,1,0);
            
            th_start = theta_bolt(1,1);
            th_ideal_step = 2*pi/(num_antennas/2); % ideal distribution is even 30deg spacing (pi/6), referenced to first bolt hole for 24 antennas, 22.5 deg for 32 antennas
            th_cables = th_start:th_ideal_step:th_start + 2*pi - th_ideal_step;
            r_cables = (tru_diam_ft/2)*ones(size(th_cables)) - 4/12; % shift in by 4"
            [x_cables, y_cables] = pol2cart(th_cables, r_cables);
            xy_vals = th_cables.*(tru_diam_ft/2);
            
            set(0,'CurrentFigure',h_2d);
            r = tru_diam_ft/2;
            conv = 180/(pi*r);
            
%             for j = 1:num_cables
%                 plot([xy_vals(j)*conv,xy_vals(j)*conv],[gbindata.eff_eave_height_ft(i), gbindata.eff_eave_height_ft(i) - gbindata.eff_eave_to_wirerope_bottom_mm(i)/304.8],'LineWidth',2,'Color',plot_cols(1,:));
%             end
            scatter(xy_vals*conv,reref_antenna_positions_ft(1:2:num_antennas-1),120,'Filled','s','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
            scatter(xy_vals*conv,reref_antenna_positions_ft(2:2:num_antennas),120,'Filled','d','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
            
            %
            %         t1 = repmat(theta_bolt*tru_diam_ft/2,length(ring_centers_ft),1);
            %         m = length(theta_bolt);
            %         t2 = repmat(ring_centers_ft,1,m);
            %         scatter(t1(:)*conv,t2(:),20,'d','r','Filled')
            %         title(rbindata.gvgsi_partno{i});
            %
            
        end
        
        ip_mtx = repmat(reref_antenna_positions_ft',1,length(ring_centers_ft)); % make into matrix
        rc_mtx = repmat(ring_centers_ft,length(reref_antenna_positions_ft'),1);
        delta = abs(ip_mtx - rc_mtx);
        [mins,idx] = min(delta,[],2); % find the closest ring center
        
        idx(idx==0) = idx(idx==0)+1;
        idx(idx>length(ring_centers_ft)) = idx(idx>length(ring_centers_ft))-1;
        
        sensor_position_ft(i,:) = ring_centers_ft(idx); % change antenna position to relevant sheet center, referenced to bottom bolt hole
        gbindata.sensor_position_mm(i,:) = (gbindata.eff_eave_height_ft(i)-sensor_position_ft(i,:))*304.8; % reref to top bolt hole, mm
    end
end

odd_idx = [1:2:23];
even_idx = [2:2:24];

% % plot all distributions
% for j = 1:num_hgroups
%     gdx = find(gbindata.height_group==j);
%     i = gdx(1);
%     
%     if (simple_plot)
%         
%         if (gbindata.corrugation(i) == 'Wide')
%             psh_ft = w_sheet_height_ft;
%             ttle = 'Wide';
%         elseif (gbindata.corrugation(i) == 'Narrow')
%             psh_ft = n_sheet_height_ft;
%             ttle = 'Narrow';
%         else
%             psh_ft = b_sheet_height_ft;
%             ttle = 'Behlen';
%         end
%         plot_cols = parula(10);
%         
%         figure;
%         hold on;
%         line([0 13],[0, 0],'Color','blue','LineWidth',1.5);
%         line([0 13],[plenum_height_ft, plenum_height_ft],'Color','blue','LineWidth',1.5);
%         line([0 13],[gbindata.eff_eave_height_ft(i),gbindata.eff_eave_height_ft(i)],'Color','blue','LineWidth',1.5);
%         for k = 1:gbindata.num_rings(i)
%             line([0 13],[psh_ft*k,psh_ft*k],'Color','black');
%         end
%         
%         scatter([1:1:12],sensor_position_ft(i,odd_idx),'Filled','s','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
%         scatter([1:1:12],sensor_position_ft(i,even_idx),'Filled','d','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
%         
%         % grid on;
%         xlim([0 13]);
%         xticks([1:1:12])
%         yticks([psh_ft/2:psh_ft:(gbindata.num_rings(i)*psh_ft - psh_ft/2)])
%         ylim([gbindata.eff_eave_height_ft(i)*-0.1 gbindata.eff_eave_height_ft(i)*1.1])
%         xlabel('Cable Number');
%         ylabel('Ring ID [ft]');
% %         title(strcat(ttle," Corrugation, ", num2str(gbindata.num_rings(i))," Rings"));
%         
%         yt = gbindata.num_rings(i):-1:1;
%         yticklabels(compose('%d',yt));
%         
%     end
%     
% end


%% evaluate positions

if (evaluate_position_metrics)
    for i = 1:length(sensor_rings.wide.nrings)
        tt = unique(sensor_rings.wide.ring_id(:,i));
        sensor_rings.wide.sens_all_rings(i) = 0;
        if (length(tt) == (sensor_rings.wide.nrings(i)-2))
            if (tt == [2:sensor_rings.wide.nrings(i)-1]')
                sensor_rings.wide.sens_all_rings(i) = 1;
            end
        end
        
        ttt = [sensor_rings.wide.ring_id(odd_idx,i),sensor_rings.wide.ring_id(even_idx,i)];
        [a,~] = size(unique(ttt,'rows'));
        if (a<12)
            sensor_rings.wide.cbl_all_unique(i) = 0;
        else
            sensor_rings.wide.cbl_all_unique(i) = 1;
        end
        
    end
    
    for i = 1:length(sensor_rings.narrow.nrings)
        tt = unique(sensor_rings.narrow.ring_id(:,i));
        sensor_rings.narrow.sens_all_rings(i) = 0;
        if (length(tt) == (sensor_rings.narrow.nrings(i)-3))
            if (tt == [2:sensor_rings.narrow.nrings(i)-2]')
                sensor_rings.narrow.sens_all_rings(i) = 1;
            end
        end
        
        ttt = [sensor_rings.narrow.ring_id(odd_idx,i),sensor_rings.narrow.ring_id(even_idx,i)];
        [a,~] = size(unique(ttt,'rows'));
        if (a<12)
            sensor_rings.narrow.cbl_all_unique(i) = 0;
        else
            sensor_rings.narrow.cbl_all_unique(i) = 1;
        end
    end
    
    for i = 1:length(sensor_rings.behlen.nrings)
        tt = unique(sensor_rings.behlen.ring_id(:,i));
        sensor_rings.behlen.sens_all_rings(i) = 0;
        if (length(tt) == (sensor_rings.behlen.nrings(i)-3))
            if (tt == [2:sensor_rings.behlen.nrings(i)-2]')
                sensor_rings.behlen.sens_all_rings(i) = 1;
            end
        end
        
        ttt = [sensor_rings.behlen.ring_id(odd_idx,i),sensor_rings.behlen.ring_id(even_idx,i)];
        [a,~] = size(unique(ttt,'rows'));
        if (a<12)
            sensor_rings.behlen.cbl_all_unique(i) = 0;
        else
            sensor_rings.behlen.cbl_all_unique(i) = 1;
        end
    end
    
    figure
    hold on;
    scatter(-10,-10,'o','Filled','b')
    scatter(-10,-10,'o','Filled','k')
    scatter(-10,-10,'o','Filled','r')
    scatter(-10,-10,'d','Filled','b')
    scatter(-10,-10,'d','Filled','k')
    scatter(-10,-10,'d','Filled','r')
    scatter(-10,-10,'s','Filled','b')
    scatter(-10,-10,'s','Filled','k')
    scatter(-10,-10,'s','Filled','r')
    
    for i = 1:length(rbindata.corrugation)
        if (rbindata.corrugation(i) == 'Wide')
            cbl_unique = sensor_rings.wide.cbl_all_unique(find(sensor_rings.wide.nrings==rbindata.num_rings(i)));
            sens_all_rings = sensor_rings.wide.sens_all_rings(find(sensor_rings.wide.nrings==rbindata.num_rings(i)));
            
            if (cbl_unique & sens_all_rings)
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'o','Filled','b')
            elseif (~cbl_unique & sens_all_rings)
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'o','Filled','k')
            else
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'o','Filled','r')
            end
        elseif (rbindata.corrugation(i) == 'Narrow')
            cbl_unique = sensor_rings.narrow.cbl_all_unique(find(sensor_rings.narrow.nrings==rbindata.num_rings(i)));
            sens_all_rings = sensor_rings.narrow.sens_all_rings(find(sensor_rings.narrow.nrings==rbindata.num_rings(i)));
            
            if (cbl_unique & sens_all_rings)
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'d','Filled','b')
            elseif (~cbl_unique & sens_all_rings)
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'d','Filled','k')
            else
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'d','Filled','r')
            end
        else
            cbl_unique = sensor_rings.behlen.cbl_all_unique(find(sensor_rings.behlen.nrings==rbindata.num_rings(i)));
            sens_all_rings = sensor_rings.behlen.sens_all_rings(find(sensor_rings.behlen.nrings==rbindata.num_rings(i)));
            
            if (cbl_unique & sens_all_rings)
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'s','Filled','b')
            elseif (~cbl_unique & sens_all_rings)
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'s','Filled','k')
            else
                scatter(rbindata.diameter_ft(i),rbindata.eave_height_ft(i),'s','Filled','r')
            end
        end
    end
    xlabel('Diameter [ft]');
    ylabel('Eave Height [ft]');
    title('Bin Sizes');
    legend('WCorr, Unique & All Rings','WCorr, Not Unique & All Rings', 'WCorr, Unique & Not All Rings','NCorr, Unique & All Rings','NCorr, Not Unique & All Rings', 'NCorr, Unique & Not All Rings','BCorr, Unique & All Rings','BCorr, Not Unique & All Rings', 'BCorr, Unique & Not All Rings','Location','SouthEast');
    xlim([0 160]);
    ylim([0 110]);
    grid on;
end


%% save outputs

if (save_structs)
    outputFolder = fullfile(pwd, 'bin_structs');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    rbindata.version = version;
    gbindata.version = version;
    save(strcat(outputFolder,'\rbindata.mat'),'rbindata');
    save(strcat(outputFolder,'\gbindata.mat'),'gbindata');
end

toc;
