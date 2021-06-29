%--------------------------------------------------------------------------

%{
CHANGES TO MAKE:

Top Box on Bottom
- Maybe modify install ground coax to be by size group instead of
  individual bin
%}


%{
Create Manufacturing, Numerological, and Installation Structs

151 Research Inc.
Kyle Nemez

Important structs built by this script:

Manufacturing Data - calculated according to bin group
mbindata
    .ferrule_positions_mm - double; row vector containing the distance from
        the top of the wire rope loop to the ferrule crimp position [mm]
    .wrope_loop_height_mm - double; height of rope stop crimped wire rope
        loop, distance is from top hanger bracket to bottom anchor bracket
        [mm]
    .wrope_cut_length_mm - double; cut length of wire rope to form a loop
        of the associated wrope_loop_height_mm [mm]
    .tbot_wall_coax_length_mm - double; cut length of coax to run from 
        sensor to top hanger bracket [mm]
    .tbob_wall_coax_length_mm - double; cut length of coax to run from 
        sensor to bottom hanger bracket [mm]
    .tbot_zzz_lengths - struct; contains sub-arrays of the manufacturing 
        lengths of the various bin size dependent non-coaxial cabling
            .cbl_002049 - double; cut length of top to bottom box cable,
                part number 002049-ZZZ-rev [mm]
            .cbl_002053 - double; cut length of top box to temp/humid cable
                extension cable, part number 002053-ZZZ-rev [mm]
            .cbl_002059 - double; cut length of the top box to peak sonar
                extension cable, part number 002059-ZZZ-rev [mm]
            .cbl_002060 - double; cut length of the top box to eave sonar
                extension cable, part number 002060-ZZZ-rev [mm]
    .tbob_zzz_lengths - struct; contains sub-arrays of the manufacturing 
        lengths of the various bin size dependent non-coaxial cabling
            .cbl_002049 - double; cut length of top to bottom box cable,
                part number 002049-ZZZ-rev [mm]
            .cbl_002053 - double; cut length of top box to temp/humid cable
                extension cable, part number 002053-ZZZ-rev [mm]

Numerical Data - calculated for each unique bin
nbindata
    .th_all_bolts - cell; each cell entry contains a row vector of the
        angular position of every bolt hole in the bin [rad]
    .tru_diam_mm - double; true diameter of the bin, not equal to the
        typically listed nominal diameter [mm]
    .th_cables - double; angular position of each assembled cable string,
        this is the estimated placement, installed location may be
        different [rad]
    .corrected_sensor_position_mm - double; distance from the true eave
        height (not bolt holes) to the center of each sensor, expressed as
        a negative number [mm]

Installation Data - calculated for each unique bin
ibindata
    .absolute_bolts - integer; nbindata.th_cables translated to absolute
        number of bolts, counted from a chosen reference bolt (bolt #1)
        clockwise when viewed from above
    .sht_bolts - integer; absolute bolts translated to bolts from the edge
        of the horizontal sheet the cable hangs in, again counted clockwise
        when viewed from above
    .sht_num - integer; number of the horizontal sheet that the cable hangs
        in, counted cfrom a chosen reference sheet (sheet #1) clockwise
        when viewed from above
    .ring_num - integer; index of the ring that each sensor resides in,
        where the uppermost ring is ring #1
    .tbot_roof_coax_lengths_mm - double; estimated cut length of each field
        terminated run of coax from the hanger bracket to the top box,
        taking the pre-defined route from bracket to top box [mm]
    .tbot_total_roof_coax_mm - double; sum of all roof_coax_lengths_mm for each
        bin, the total amount of coax that must be shipped to an install
        site [mm]
    .tbob_roof_coax_lengths_mm - double; estimated cut length of each field
        terminated run of coax from the hanger bracket to the top box,
        taking the pre-defined route from bracket to top box [mm]
%}
%--------------------------------------------------------------------------

tic;

clear all;
close all;
clc;

version = '1.2.6';

plot_single_array = 0;
calculate_tbobvstbot = 0;
save_structs = 1;

base_path = pwd; % get file path to script

load(strcat(base_path,'\bin_structs\rbindata.mat'));
load(strcat(base_path,'\bin_structs\gbindata.mat'));
addpath(strcat(base_path,'\sub_functions'));

%% calculate manufacturing outputs

%--------------------------------------------------------------------------
% ANY SIGNFICANT EDITS IN THIS SECTION SHOULD RESULT IN AN UPDATE TO THE
% VERSION NUMBER CONTAINED IN 'version'
%--------------------------------------------------------------------------

% important build parameters

% bin parameters
floor_clearance_ft = 0.762/0.3048; % clearance from plenum floor to lowest  bottom of sensor
plenum_height_ft = 18/12; %UPDATE 6/10 KD - changed to 18" worst case plenum flor height is 17.25", this is referenced to lowest ring of bolt holes on bottom sheet

% distance from eave bolt ring to bottom of hole wire rope is pulled
% through
% columns represent different corrugations:
% (1,:) -> Wide corrugation (4")
% (2,:) -> Narrow corrugation (2.66")
% (3,:) -> Behlen corruagation (3.125")
% rows represent mounting on the highest available corruagtion or the next
% one down
% (:,1) -> highest available peak
% (:,2) -> 2nd highest available peak
bolt_to_brackethole_mm = zeros(3,2);
bolt_to_brackethole_mm(1,:) = [115.86,217.49];
bolt_to_brackethole_mm(2,:) = [81.854,149.418];
bolt_to_brackethole_mm(3,:) = [93.665,173.04];

brackethole_to_hanger_mm = 49.73; % distance from bottom of hole to the hang point on the bracket for the wire rope [mm]

overlap_mm = 20; % overlap for wire rope loop crimping
wrfudge_mm = 20; % fudge factor for accounting for loop ends/curvature

% cumulative shift for antennas wrt eave bolts is down
% bolt_to_brackethole_mm and up brackethole_to_hanger_mm
% assume that the install is on highest peak available, 2nd highest peak
% data is just for reference when an install varies
net_shift_mm = zeros(size(gbindata.antenna_length_ft));

wg_idx = find(gbindata.corrugation == 'Wide');
ng_idx = find(gbindata.corrugation == 'Narrow');
bg_idx = find(gbindata.corrugation == 'Behlen');

net_shift_mm(wg_idx) = bolt_to_brackethole_mm(1,1) - brackethole_to_hanger_mm;
net_shift_mm(ng_idx) = bolt_to_brackethole_mm(2,1) - brackethole_to_hanger_mm;
net_shift_mm(bg_idx) = bolt_to_brackethole_mm(3,1) - brackethole_to_hanger_mm;

% ferrule positions are sensor position minus bolt to brackethole, minus
% half an antenna length, plus the brackethole to hanger distance
% ferrule position is referenced to the peak of the wire rope loop
mbindata.ferrule_positions_mm = gbindata.sensor_position_mm - (gbindata.antenna_length_ft*304.8/2) - net_shift_mm;
% mbindata.ferrule_positions_mm = gbindata.sensor_position_mm - (bolt_to_brackethole_mm + (gbindata.antenna_length_ft*304.8/2)) + brackethole_to_hanger_mm;

% wire rope loop height is eff eave height to bottom minus bolt to
% brackethole, plus brackethole to hanger, rounded to nearest 10mm. 
gbindata.eff_eave_to_wirerope_bottom_mm = ((gbindata.eff_eave_height_ft - plenum_height_ft) - floor_clearance_ft + ((7+6)/12))*304.8;
mbindata.wrope_loop_height_mm = 10*ceil((gbindata.eff_eave_to_wirerope_bottom_mm - net_shift_mm)/10);
% mbindata.wrope_loop_height_mm = 10*ceil((gbindata.eff_eave_to_wirerope_bottom_mm - bolt_to_brackethole_mm + brackethole_to_hanger_mm)/10);

% cut length is twice the height plus overlap plus fudge
mbindata.wrope_cut_length_mm = 2*mbindata.wrope_loop_height_mm + overlap_mm + wrfudge_mm;

% determine length of coax from antenna to mounting bracket
wcoax_extra_mm = 300; % add 300mm of extra coax

% TBoT coax from antenna to bracket (wall coax) is ferrule positions plus extra
mbindata.tbot_wall_coax_length_mm = mbindata.ferrule_positions_mm + wcoax_extra_mm;
mbindata.tbot_wall_coax_length_mm = 10*ceil(mbindata.tbot_wall_coax_length_mm/10); % round up to nearest cm

% TBoB coax length is sensor position minus 1/2 antenna length to wire rope
% bottom plus extra
mbindata.tbob_wall_coax_length_mm = repmat(gbindata.eff_eave_height_ft*304.8,1,24) - (gbindata.sensor_position_mm + (gbindata.antenna_length_ft*304.8/2));
% mbindata.tbob_wall_coax_length_mm = repmat(gbindata.eff_eave_to_wirerope_bottom_mm,1,24) + wcoax_extra_mm - (gbindata.sensor_position_mm + (gbindata.antenna_length_ft*304.8/2));
mbindata.tbob_wall_coax_length_mm = 10*ceil(mbindata.tbob_wall_coax_length_mm/10); % round up to nearest cm

% calculate zzz lengths for 6 conductor (TBoT install)
roof_fudge = 1; % multiplier for roof length
peak_to_tbox_ft = 12; % distance from bin peak to top box install location
mbindata.tbot_zzz_lengths = calculate_tbot_zzz_length(gbindata.eff_eave_height_ft,gbindata.max_gdiam_ft,30,roof_fudge,peak_to_tbox_ft);

% calculate zzz lengths for 6 conductor (TBoB install)
tbob_roof_fudge = 1; % multiplier for roof length
tbox_fudge_length_ft = 4; % additive (or subtractive) parameter to adjust length of t/rh extension cable
mbindata.tbob_zzz_lengths = calculate_tbob_zzz_length(gbindata.eff_eave_height_ft,gbindata.max_gdiam_ft,30,tbob_roof_fudge,tbox_fudge_length_ft);

%% calculate installation and numerical outputs

%--------------------------------------------------------------------------
% ANY SIGNFICANT EDITS IN THIS SECTION SHOULD RESULT IN AN UPDATE TO THE
% VERSION NUMBER CONTAINED IN 'version'
%--------------------------------------------------------------------------

% this section is unique for every bin.
[~,num_antennas] = size(mbindata.ferrule_positions_mm);
num_cables = num_antennas/2;
num_bins = numel(rbindata.gvgsi_partno);
nbindata.th_all_bolts = cell(num_bins,1);
nbindata.tru_diam_mm = zeros(num_bins,1);
nbindata.th_cables = zeros(num_bins,num_antennas/2);

% bin configuration parameters set for worst case
stiff_per_sheet = 3;
alt_stagger = 1;

plot_2d = 0;
plot_3d = 0;
plot_cols = parula(10);

% obtain all bolt angles and assign cables to angular bolt positions
for i = 1:num_bins

    % get parameters of the bin
    [h_2d,h_3d,tru_diam_ft,sheet_centers_ft,theta_bolt,theta_all_bolts] = plot_bin_panels(rbindata.num_rings(i),rbindata.corrugation(i),rbindata.diameter_ft(i),plenum_height_ft,stiff_per_sheet,alt_stagger,plot_2d,plot_3d);
    
    nbindata.th_all_bolts{i} = theta_all_bolts; 
    nbindata.tru_diam_mm(i,1) = tru_diam_ft*304.8;
    
    % adjust angles to be on bolts
    th_start = theta_bolt(1,1);
    th_ideal_step = 2*pi/(num_antennas/2); % ideal distribution is even 30deg spacing (pi/6), referenced to first bolt hole for 24 antennas, 22.5 deg for 32 antennas
    th_ideal = th_start:th_ideal_step:th_start + 2*pi - th_ideal_step; 
    th_ideal_mtx = repmat(th_ideal',1,length(theta_bolt)); % make into matrix
    th_bolt_mtx = repmat(theta_bolt,length(th_ideal'),1);
    th_delta = abs(th_ideal_mtx - th_bolt_mtx);
    [th_mins,th_idx] = min(th_delta,[],2); % find the closest bolt hole

    nbindata.th_cables(i,:) = theta_bolt(th_idx);
    % add half of the difference between actual height and effective height
    % to the grouped sensor position to account for the distance between
    % the top bolt holes and the true top of the bin
    % also make it a negative number
    nbindata.corrected_sensor_position_mm(i,:) = -1*(gbindata.sensor_position_mm(rbindata.size_code(i),:) + (rbindata.eave_height_ft(i) - gbindata.eff_eave_height_ft(rbindata.size_code(i)))*304.8/2);
    
    % plot to verify
    if (h_2d ~= 0) || (h_3d ~= 0)
        th_cables = nbindata.th_cables(i,:);
        r_cables = (tru_diam_ft/2)*ones(size(th_cables)) - 4/12; % shift in by 4"
        [x_cables, y_cables] = pol2cart(th_cables, r_cables);
        xy_vals = th_cables.*(tru_diam_ft/2);
        if (h_2d ~= 0)
            set(0,'CurrentFigure',h_2d);
            r = tru_diam_ft/2;
            conv = 180/(pi*r);

            for j = 1:num_cables
                plot([xy_vals(j)*conv,xy_vals(j)*conv],[gbindata.eff_eave_height_ft(rbindata.size_code(i)), gbindata.eff_eave_height_ft(rbindata.size_code(i)) - gbindata.eff_eave_to_wirerope_bottom_mm(rbindata.size_code(i))/304.8],'LineWidth',2,'Color',plot_cols(1,:));
            end
            scatter(xy_vals*conv,gbindata.eff_eave_height_ft(rbindata.size_code(i)) - gbindata.sensor_position_mm(rbindata.size_code(i),1:2:num_antennas-1)/304.8,120,'Filled','s','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
            scatter(xy_vals*conv,gbindata.eff_eave_height_ft(rbindata.size_code(i)) - gbindata.sensor_position_mm(rbindata.size_code(i),2:2:num_antennas)/304.8,120,'Filled','d','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);


            t1 = repmat(theta_bolt*tru_diam_ft/2,length(sheet_centers_ft),1);
            m = length(theta_bolt);
            t2 = repmat(sheet_centers_ft,1,m);
            scatter(t1(:)*conv,t2(:),20,'d','r','Filled')
            title(rbindata.gvgsi_partno{i});
        end
        if (h_3d ~= 0)
            set(0,'CurrentFigure',h_3d);

            for j = 1:num_cables
                plot3([x_cables(j), x_cables(j)],[y_cables(j), y_cables(j)],[gbindata.eff_eave_height_ft(rbindata.size_code(i)), gbindata.eff_eave_height_ft(rbindata.size_code(i)) - gbindata.eff_eave_to_wirerope_bottom_mm(rbindata.size_code(i))/304.8],'LineWidth',2,'Color',plot_cols(1,:));
            end
            scatter3(x_cables', y_cables',gbindata.eff_eave_height_ft(rbindata.size_code(i)) - gbindata.sensor_position_mm(rbindata.size_code(i),1:2:num_antennas-1)/304.8,120,'Filled','s','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);
            scatter3(x_cables', y_cables',gbindata.eff_eave_height_ft(rbindata.size_code(i)) - gbindata.sensor_position_mm(rbindata.size_code(i),2:2:num_antennas)/304.8,120,'Filled','d','MarkerFaceColor',plot_cols(10,:),'MarkerEdgeColor',plot_cols(1,:),'LineWidth',2);

            [x_bolts, y_bolts] = pol2cart(theta_bolt, tru_diam_ft/2*ones(size(theta_bolt)));
            x1 = repmat(x_bolts,length(sheet_centers_ft),1);
            y1 = repmat(y_bolts,length(sheet_centers_ft),1);
            m = length(x_bolts);
            t3 = repmat(sheet_centers_ft,1,m);
            scatter3(x1(:),y1(:),t3(:),20,'d','r','Filled');
            title(rbindata.gvgsi_partno{i});
        end
    end

end

% get install parameters
bolt_holes_per_sheet = 12*ones(size(nbindata.tru_diam_mm));
bolt_holes_per_sheet(rbindata.corrugation=='Behlen') = 32;
bolt_spacing = (9+3/8)*ones(size(nbindata.tru_diam_mm));
bolt_spacing(rbindata.corrugation=='Behlen') = (125.984)/32;

% translate angular position to actual bolt count
ibindata.absolute_bolts = round((nbindata.th_cables(:,1:num_cables).*nbindata.tru_diam_mm/(2*304.8)).*(12./(bolt_spacing)),0);
ibindata.sht_bolts = mod(ibindata.absolute_bolts,bolt_holes_per_sheet);
ibindata.sht_num = floor(ibindata.absolute_bolts./bolt_holes_per_sheet) + 1;

% translate distance from eave bolt holes to ring number
wide_idx = find(rbindata.corrugation == 'Wide');
narrow_idx = find(rbindata.corrugation == 'Narrow');
behlen_idx = find(rbindata.corrugation == 'Behlen');
for i = 1:length(wide_idx)
    j = wide_idx(i);
    ibindata.ring_num(j,:) = round(2*(0.5 + (gbindata.sensor_position_mm(rbindata.size_code(j),:)/304.8)/(44/12)),0)/2;
end
for i = 1:length(narrow_idx)
    j = narrow_idx(i);
    ibindata.ring_num(j,:) = round(0.5 + (gbindata.sensor_position_mm(rbindata.size_code(j),:)/304.8)/(32/12),0);
end
for i = 1:length(behlen_idx)
    j = behlen_idx(i);
    ibindata.ring_num(j,:) = round(0.5 + (gbindata.sensor_position_mm(rbindata.size_code(j),:)/304.8)/(39.3701/12),0);
end

% calculate cut lengths and total length required for field terminated
% roof coax
cable_run_diameter_ft = 5; % set diameter of cable run around roof cap
ibindata.tbot_roof_coax_lengths_mm = calculate_roof_coax(rbindata.diameter_ft,rbindata.roof_angle,roof_fudge,peak_to_tbox_ft,cable_run_diameter_ft,nbindata.th_cables); 
ibindata.tbot_roof_coax_lengths_mm = 500*ceil(ibindata.tbot_roof_coax_lengths_mm/500); % round up to nearest 500mm
ibindata.tbot_total_roof_coax_mm = sum(ibindata.tbot_roof_coax_lengths_mm,2);

tb_height_mm = 1500; % set height of top box above ground
for i = 1:num_bins
    tb_angle = nbindata.th_all_bolts{i}(1+bolt_holes_per_sheet(i))/2; % angle is 1/2 way around the first sheet
%     tb_angle = nbindata.th_cables(i,1); % at cable 1 for now.
    ibindata.tbob_ground_coax_length_mm(i,:) = calculate_ground_coax(rbindata.diameter_ft(i), nbindata.th_cables(i,:), rbindata.corrugation(i), tb_angle, tb_height_mm);
end

%% save output structs
if (save_structs)
    outputFolder = fullfile(pwd, 'mni_structs');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    rbindata.version = version;
    gbindata.version = version;
    save(strcat(outputFolder,"\mbindata.mat"),'mbindata', 'net_shift_mm', 'bolt_to_brackethole_mm','brackethole_to_hanger_mm');
    save(strcat(outputFolder,"\nbindata.mat"),'nbindata');
    save(strcat(outputFolder,"\ibindata.mat"),'ibindata');
end

%% calc coax attenuation, in db/mm

if (calculate_tbobvstbot)
    f = [1:1:1000];
    coax.raw.deca60030(1,:) = [10, 100, 200, 400, 1000];
    coax.raw.deca60030(2,:) = [1.1,3.8,5.6,8.4,14.5]/30480;
    p = polyfit(coax.raw.deca60030(1,:),coax.raw.deca60030(2,:),3);
    coax.deca60030 = zeros(size(f));
    for i = 1:1:3+1
        coax.deca60030 = coax.deca60030 + p(i)*f.^(3+1-i);  % pretty good curve fit below 500MHz
    end
    coax.deca60030 = -1*coax.deca60030;
    plot_data = true;
    if (plot_data)
        figure
        plot(coax.raw.deca60030(1,:),-1*coax.raw.deca60030(2,:));
        hold on;
        plot(f,coax.deca60030);
        %     xlim([10 100]);
        grid on;
        xlabel('Frequency [MHz]');
        ylabel('Cable Attenuation [dB/mm]');
        legend('Datasheet','Curve');
        title('DECA 60-030');
    end
    
    deca_pul_losses(1,:) = [10,30,100,500,1000];
    deca_pul_losses(2,:) = coax.deca60030(1,deca_pul_losses(1,:));
    
    
    odd_idx = 1:2:num_antennas-1;
    even_idx = 2:2:num_antennas;
    for i = 1:num_bins
        coax_length_mm(i,:) = ibindata.tbot_roof_coax_lengths_mm(i,:) + mbindata.tbot_wall_coax_length_mm(rbindata.size_code(i),:);
        tbob_coax_length_mm(i,:) = ibindata.tbob_ground_coax_length_mm(i,:) + mbindata.tbob_wall_coax_length_mm(rbindata.size_code(i),:);
        delta_coax_length_mm(i,odd_idx) = tbob_coax_length_mm(i,odd_idx) - coax_length_mm(i,even_idx);
        delta_coax_length_mm(i,even_idx) = tbob_coax_length_mm(i,even_idx) - coax_length_mm(i,odd_idx);
        delta_nwall_coax_length_mm(i,:) = ibindata.tbob_ground_coax_length_mm(i,:) - ibindata.tbot_roof_coax_lengths_mm(i,:);
        delta_losses(i,:,:) = (deca_pul_losses(2,:).'*delta_coax_length_mm(i,:)).';
        delta_nwall_losses(i,:,:) = (deca_pul_losses(2,:).'*(ibindata.tbob_ground_coax_length_mm(i,:) - ibindata.tbot_roof_coax_lengths_mm(i,:))).';
    end
    
    vol = (pi*(rbindata.diameter_ft*0.3048/2).^2).*rbindata.eave_height_ft*0.3048;
    [~,I] = sort(vol);
    
    num_groups = length(I);
    hcol = jet(num_groups);
    
    alpha_name = {'Gilmore','Semler','Binscarth','Alix','Record','Agtegra','I-80','CVA','Balgonie','Azteca'};
    alpha_diam_ft = [32.8,47.9,78,53.7,60,135,105,105,78,90];
    alpha_height_ft = [32,62.6,84.5,44.2,51.5,88,74.8,98.5,106.5,61.3];
    alpha_vol = pi*((0.3048*alpha_diam_ft/2).^2).*(0.3048*alpha_height_ft);
    for i = 1:length(alpha_vol)
        [k(i), j(i)] = min(abs(alpha_vol(i)-vol));
    end
    
    
    figure
    subplot(1,3,1)
    hold on;
    for i = 1:num_groups
        plot(1:12,ibindata.tbot_roof_coax_lengths_mm(I(i),odd_idx)'/1000,'Color',hcol(i,:));
    end
    xlim([1 12]);
    ylim([1 80]);
    grid on;
    xlabel('Cable Number');
    ylabel('Coax Length [m]');
    title('TBoT Roof Coax Length for All Bins');
    
    subplot(1,3,2)
    hold on;
    for i = 1:num_groups
        plot(1:12,ibindata.tbob_ground_coax_length_mm(I(i),odd_idx)'/1000,'Color',hcol(i,:));
    end
    xlim([1 12]);
    ylim([1 80]);
    grid on;
    xlabel('Cable Number');
    ylabel('Coax Length [m]');
    title('TBoB Ground Coax Length for All Bins');
    
    subplot(1,3,3)
    hold on;
    for i = 1:num_groups
        plot(1:12,delta_nwall_coax_length_mm(I(i),odd_idx)'/1000,'Color',hcol(i,:));
    end
    xlim([1 12]);
    grid on;
    xlabel('Cable Number');
    ylabel('Delta Coax Length [m]');
    title('(TBoB-TBoT) Coax Length for All Bins');
    
    
    
    figure
    subplot(1,3,1)
    hold on;
    for i = 1:num_groups
        plot(1:24,coax_length_mm(I(i),:)'/1000,'Color',hcol(i,:));
    end
    xlim([1 24]);
    ylim([1 100]);
    grid on;
    xlabel('Antenna Number');
    ylabel('Coax Length [m]');
    title('TBoT Coax Length for All Bins');
    
    subplot(1,3,2)
    hold on;
    for i = 1:num_groups
        plot(1:24,tbob_coax_length_mm(I(i),:)'/1000,'Color',hcol(i,:));
    end
    xlim([1 24]);
    ylim([1 100]);
    grid on;
    xlabel('Antenna Number');
    ylabel('Coax Length [m]');
    title('TBoB Coax Length for All Bins');
    
    subplot(1,3,3)
    hold on;
    for i = 1:num_groups
        plot(1:24,delta_coax_length_mm(I(i),:)'/1000,'Color',hcol(i,:));
    end
    xlim([1 24]);
    grid on;
    xlabel('Antenna Number');
    ylabel('Delta Coax Length [m]');
    title('(TBoB-TBoT) Coax Length for All Bins');
    
    figure
    for i = 1:length(deca_pul_losses(1,:))
        subplot(1,length(deca_pul_losses(1,:)),i);
        plot(1:24,squeeze(delta_losses(:,:,i)))
        xlim([1 24]);
        ylim([-20 20]);
        grid on;
        xlabel('Antenna Number');
        ylabel('Signal Level Change [dB]');
        title(strcat('Signal Level Relative To TBoT, f = ',num2str(deca_pul_losses(1,i)),' [MHz]'));
    end
    
    figure
    subplot(1,2,1);
    scatter(vol,sum(delta_coax_length_mm/1000,2),'MarkerEdgeColor','c','MarkerFaceColor','c')
    hold on;
    for i = 1:length(alpha_vol)
        scatter(vol(j(i)),sum(delta_coax_length_mm(j(i),:))/1000,'k','Filled');
        text(vol(j(i))+50,3+sum(delta_coax_length_mm(j(i),:))/1000,alpha_name{i});
    end
    xlabel('Bin Volume [m^3]');
    ylabel('Differential Total Coax [m]');
    title('Relative Total Coax Usage for TBoB Compared to TBoT');
    set(gca,'xscale','log')
    grid on;
    
    subplot(1,2,2);
    hold on;
    scatter(vol(1),mean(delta_coax_length_mm(1,:)/1000,2),'MarkerEdgeColor','c','MarkerFaceColor','c');
    scatter(vol(1),max(delta_coax_length_mm(1,:)/1000,[],2),'MarkerEdgeColor','m','MarkerFaceColor','m');
    scatter(vol(1),min(delta_coax_length_mm(1,:)/1000,[],2),'MarkerEdgeColor','g','MarkerFaceColor','g');
    
    scatter(vol,mean(delta_coax_length_mm/1000,2),'MarkerEdgeColor','c','MarkerFaceColor','c');
    scatter(vol,max(delta_coax_length_mm/1000,[],2),'MarkerEdgeColor','m','MarkerFaceColor','m');
    scatter(vol,min(delta_coax_length_mm/1000,[],2),'MarkerEdgeColor','g','MarkerFaceColor','g');
    
    for i = 1:length(alpha_vol)
        line([vol(j(i)) vol(j(i))],[min(delta_coax_length_mm(j(i),:))/1000 max(delta_coax_length_mm(j(i),:))/1000],'Color',[0.5 0.5 0.5])
        scatter(vol(j(i)),mean(delta_coax_length_mm(j(i),:))/1000,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(vol(j(i)),max(delta_coax_length_mm(j(i),:))/1000,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(vol(j(i)),min(delta_coax_length_mm(j(i),:))/1000,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        text(vol(j(i))+20,1+mean(delta_coax_length_mm(j(i),:))/1000,alpha_name{i});
    end
    ylim([-45 45])
    xlabel('Bin Volume [m^3]');
    ylabel('Coax Length Change [m]');
    title('Relative Coax Length for TBoB Compared to TBoT');
    set(gca,'xscale','log')
    grid on;
    legend('Mean Coax Length Change','Max Coax Length Change','Min Coax Length Change','Location','NorthWest');
    
    
    
    
    
    delta_coax_length_mm2 = tbob_coax_length_mm - coax_length_mm;
    figure
    subplot(1,2,1);
    scatter(vol,sum(delta_coax_length_mm2/1000,2),'MarkerEdgeColor','c','MarkerFaceColor','c')
    hold on;
    for i = 1:length(alpha_vol)
        scatter(vol(j(i)),sum(delta_coax_length_mm2(j(i),:))/1000,'k','Filled');
        text(vol(j(i))+50,3+sum(delta_coax_length_mm2(j(i),:))/1000,alpha_name{i});
    end
    xlabel('Bin Volume [m^3]');
    ylabel('Differential Total Coax [m]');
    title('Relative Total Coax Usage for TBoB Compared to TBoT');
    set(gca,'xscale','log')
    grid on;
    
    subplot(1,2,2);
    hold on;
    scatter(vol(1),mean(delta_coax_length_mm2(1,:)/1000,2),'MarkerEdgeColor','c','MarkerFaceColor','c');
    scatter(vol(1),max(delta_coax_length_mm2(1,:)/1000,[],2),'MarkerEdgeColor','m','MarkerFaceColor','m');
    scatter(vol(1),min(delta_coax_length_mm2(1,:)/1000,[],2),'MarkerEdgeColor','g','MarkerFaceColor','g');
    
    scatter(vol,mean(delta_coax_length_mm2/1000,2),'MarkerEdgeColor','c','MarkerFaceColor','c');
    scatter(vol,max(delta_coax_length_mm2/1000,[],2),'MarkerEdgeColor','m','MarkerFaceColor','m');
    scatter(vol,min(delta_coax_length_mm2/1000,[],2),'MarkerEdgeColor','g','MarkerFaceColor','g');
    
    for i = 1:length(alpha_vol)
        line([vol(j(i)) vol(j(i))],[min(delta_coax_length_mm2(j(i),:))/1000 max(delta_coax_length_mm2(j(i),:))/1000],'Color',[0.5 0.5 0.5])
        scatter(vol(j(i)),mean(delta_coax_length_mm2(j(i),:))/1000,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(vol(j(i)),max(delta_coax_length_mm2(j(i),:))/1000,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(vol(j(i)),min(delta_coax_length_mm2(j(i),:))/1000,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        text(vol(j(i))+20,1+mean(delta_coax_length_mm2(j(i),:))/1000,alpha_name{i});
    end
    ylim([-45 45])
    xlabel('Bin Volume [m^3]');
    ylabel('Mean Additional Coax [m]');
    title('Relative Coax Length for TBoB Compared to TBoT');
    set(gca,'xscale','log')
    grid on;
    legend('Mean Coax Length Change','Max Coax Length Change','Min Coax Length Change','Location','NorthWest');
end

toc;