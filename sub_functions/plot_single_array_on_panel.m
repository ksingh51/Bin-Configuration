function [h_2d, diameter_ft, sheet_centers_ft, theta_bolt] = plot_single_array_on_panel(num_rings, corrugation, nominal_diam_ft, plenum_height_ft, stiff_per_sheet, alt_stagger, plot2d)
%{

plots the bin sheets for the input bin

inputs:
num_rings - number of vertical sheets in the bin
corrugation - sheet corrugation, 'W' or 'N'
nominal_diam_ft - nominal diameter of the bin [ft]
plenum_height_ft - distance from base of bin to plenum floor [ft]
stiff_per_sheet - number of stiffeners per sheet, 0 2 or 3
plot2d - 0 to turn off 3d plot, 1 to turn on 2d plot

outputs:
h_2d - handle for 2d figure, or zero if plot2d is zero
diameter_ft - the true diameter of the bin in feet
sheet_centers_ft - the height of the centers of each ring, referenced to
    bottom ring of bolts
theta_bolt - angles of viable bolt holes for antenna mounting

%}

% effective arc length of sheet is 112.5" (center to center of horizontally
% adjacent sheets), except Behlen bins, which is
if (strcmp('Behlen',char(corrugation)))
    sheet_arc_eff = 125.984/12; % (3.2m)
else
    sheet_arc_eff = 112.5/12; 
end

% bolt pattern is 9 3/8", except Behlen bins, which is 3-15/16" or 0.1m
if (strcmp('Behlen',char(corrugation)))
    bolt_space = sheet_arc_eff/32; 
else
    bolt_space = (9 + 3/8)/12; 
end

% determine effective sheet height, 44" for wide corr, 32" for narrow corr
% (center to center of vertically adjacent sheets)
% default is narrow corr
if (strcmp('Wide',char(corrugation)))
    sheet_height_eff = 44/12;
elseif (strcmp('Narrow',char(corrugation)))
    sheet_height_eff = 32/12;
else
    sheet_height_eff = 39.3701/12; % 1m
end

% get center of each ring, referenced to bottom ring of bolts
sheet_centers_ft = [sheet_height_eff/2:sheet_height_eff:num_rings*sheet_height_eff];

% number of horizontal sheets in a ring is nominal diameter/3
if (strcmp('Behlen',char(corrugation)))
    % for behlen it's metric
    num_sheets = round((nominal_diam_ft*0.3048),0);
else
    num_sheets = floor(nominal_diam_ft/3);
end

% true diameter = circumference / pi, circumference is num_sheets *
% effective arc length
diameter_ft = num_sheets*sheet_arc_eff/pi;

% get theta (in radians) of all possible bolt holes
delta_theta_bolt = bolt_space/(diameter_ft/2);
theta_bolt = 0:delta_theta_bolt:2*pi-delta_theta_bolt;

% get theta (in radians) of all obstacles
delta_theta_sheet = (sheet_arc_eff/(diameter_ft/2))/3;
theta_sheet = 0:delta_theta_sheet:2*pi-delta_theta_sheet; % sheet edges
theta_stiff3 = delta_theta_sheet/2:delta_theta_sheet:2*pi-delta_theta_sheet/2; % stiffeners for 3 per sheet
delta_theta_stiff2 = (sheet_arc_eff/(diameter_ft/2))/2;
theta_stiff2 = delta_theta_stiff2/2:delta_theta_stiff2:2*pi-delta_theta_stiff2/2; % stiffeners for 2 per sheet

[theta_sheet_idx,~] = ismember(round(theta_bolt,2),round(theta_sheet,2));
[theta_stiff3_idx,~] = ismember(round(theta_bolt,2),round(theta_stiff3,2));
[theta_stiff2_idx,~] = ismember(round(theta_bolt,2),round(theta_stiff2,2));
theta_bolt = theta_bolt(~(theta_sheet_idx | theta_stiff3_idx | theta_stiff2_idx));


eave_height_ft = num_rings * sheet_height_eff;

height_seed = 0;
xy_seed = 0;
xy_start = xy_seed;
edge_color = [160 160 160]./255;%[128 128 128]./255;
plenum_color = hsv(1);
stiff_color = [205,92,92]./255;%[80,80,80]./255;
edge_width = 1;
stiff_width = 1.3;
num_stiffener = stiff_per_sheet * num_sheets; % total num stiffeners
dth_stiff = 2*pi/num_stiffener; % even delta theta between them

if (plot2d & (length(num_rings) == 1))
    % generate 2d unwrapped bin plot
    h_2d = figure;
    hold on;
    for i = 1:num_rings
        % theta = s/r
        r = diameter_ft/2;
        conv = 180/(pi*r);
        % plot first vertical edge of first sheet
        plot([xy_seed*conv, xy_seed*conv],[height_seed + sheet_height_eff*(i-1), height_seed + sheet_height_eff*(i)],'Color',edge_color,'LineWidth',edge_width);
        for j = 2:num_sheets+1
            % plot bottom edge of sheet
            plot([(xy_seed + sheet_arc_eff*(j-2))*conv, (xy_seed + sheet_arc_eff*(j-1))*conv],[height_seed + sheet_height_eff*(i-1), height_seed + sheet_height_eff*(i-1)],'Color',edge_color,'LineWidth',edge_width);
            if ((i ==1) & (plenum_height_ft ~= 0))
                % plot plenum floor
                %                     plot([(xy_seed + sheet_arc_eff*(j-2))*conv, (xy_seed + sheet_arc_eff*(j-1))*conv],[height_seed + plenum_height_ft, height_seed + plenum_height_ft],'Color',plenum_color,'LineWidth',edge_width);
            end
            % plot second vertical edge of sheet
            plot([(xy_seed + sheet_arc_eff*(j-1))*conv, (xy_seed + sheet_arc_eff*(j-1))*conv],[height_seed + sheet_height_eff*(i-1), height_seed + sheet_height_eff*(i)],'Color',edge_color,'LineWidth',edge_width);
            % plot top edge of sheet
            plot([(xy_seed + sheet_arc_eff*(j-2))*conv, (xy_seed + sheet_arc_eff*(j-1))*conv],[height_seed + sheet_height_eff*(i), height_seed + sheet_height_eff*(i)],'Color',edge_color,'LineWidth',edge_width);
        end
        
        % change the xy_seed to stagger the panels
        if (stiff_per_sheet == 3)
            if (alt_stagger == 0) % 2/3 wrap
                switch xy_seed
                    case 0
                        xy_seed = sheet_arc_eff/3;
                    case sheet_arc_eff/3
                        xy_seed = -1*sheet_arc_eff/3;
                    case -1*sheet_arc_eff/3
                        xy_seed = 0;
                end
            else % switch back and forth on 1/3 wrap
                switch xy_seed
                    case 0
                        xy_seed = sheet_arc_eff/3;
                    case sheet_arc_eff/3
                        xy_seed = 0;
                end
            end
        else
            switch xy_seed
                case 0
                    xy_seed = sheet_arc_eff/2;
                case sheet_arc_eff/2
                    xy_seed = 0;
            end
        end
    end
    
    axis equal;
    if (stiff_per_sheet == 3)
        if (alt_stagger == 0)
            xlim([(-sheet_arc_eff/3 -1)*conv, (num_sheets*sheet_arc_eff + sheet_arc_eff/3 + 1)*conv]);
        else
            xlim([-1*conv, (num_sheets*sheet_arc_eff + sheet_arc_eff/3 + 1)*conv]);
        end
    else
        xlim([-1*conv, (num_sheets*sheet_arc_eff + sheet_arc_eff/2 + 1)*conv]);
    end
    ylim([-1 (eave_height_ft + 1)]);
    
    
    xlabel('Angle [degrees]');
    ylabel('Height [ft]');
    
else
    h_2d = 0;
end

end

