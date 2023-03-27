close all
clear


%% parameters
SLMm = 1152;
SLMn = 960; % half width of 1920

SLMm_mm = 10.7;
SLMn_mm = 8.8; % half width 

add_scope_params = 1;
microscope_NA_eff = [0.61, 0.565, 0.48, 0.415];
microscope_SLM_pixels = [1152, 1152, 1152, 1152];
microscope_FOV_size = [511, 511, 637.4, 637.4]*1e-6;
microscope_names = {'25X imaging', '25X stim', '20X imaging', '20X stim'};
microscope_symbols = {'*', '*', 'o', 'o'};
microscope_colors = {'g', 'r', 'g', 'r'};

obj_mag = 25;
obj_NA = 1.05;
% 
% obj_mag = 20;
% obj_NA = .95;

wavelength = 940e-9;

k = 2*pi/wavelength;

obj_RI = 1.33;
f_tube = 180;  % in mm

% efficiency parameters
max_dist = 500;
NA_eff_plot = 0.1:0.1:1;

z_range = (0:10:max_dist)*1e-6;
x_range = (0:10:max_dist)*1e-6;


NA_eff_comp = 0.1:0.01:1.05;
pat_pixels = [512, 960, 1024, 1152, 1920];

eff_lim = .95;
d_z = 2e-6;
d_x = 2e-6;


%% diffraction efficiency of blazed gratings vs number of levels

N = 1:100;% number of levels in blazed

nu = (sin(pi./N)./(pi./N)).^2;

% figure;
% subplot(1,2,1);
% plot(N, nu, 'k', 'linewidth', 2); 
% axis tight;
% title('Blazed grating diffraction efficiency')
% ylabel('diffraction efficiency');
% xlabel('Number of grating levels per period');
% 
% subplot(1,2,2);
% semilogx(N, nu, 'k', 'linewidth', 2); 
% axis tight;
% title('Blazed grating diffraction efficiency')
% ylabel('diffraction efficiency');
% xlabel('Number of grating levels per period');

figure;
semilogx(N, nu, 'k', 'linewidth', 2); 
axis tight;
title('Blazed grating diffraction efficiency')
ylabel('diffraction efficiency');
xlabel('Number of grating levels per period');


%%

M_phase = 56/55; % magnification of SLM phase to back of objective
M_image = 33/1400; % magnification of points image pattern from right after SLM

%SLM_pix_pitch = 9.2e-6;


SLM_pix_pitch = (5:0.2:12)*1e-6;

min_pix_num = [2, 4, 8, 16];
SLM_def_dist_all = zeros(numel(SLM_pix_pitch), numel(min_pix_num));

num_pix = numel(min_pix_num);

colors1 = parula(num_pix+1);

figure; hold on;
for n_px = 1:numel(min_pix_num)

    % period*sin(theta) = m*wavelength
    SLM_def_angle = asin(wavelength./(SLM_pix_pitch*min_pix_num(n_px)));

    f_lens = .3;

    SLM_def_dist = f_lens*tan(SLM_def_angle) * M_image;
    
    SLM_def_dist_all(:,n_px) = SLM_def_dist;
    
    plot(SLM_pix_pitch/1e-6, SLM_def_dist./1e-6, 'linewidth', 2, 'color', colors1(n_px,:));
end


xlabel('pixel pitch (um)')
ylabel('deflection dist (um)');
title('Blazed period vs deflection distance')
legend({[num2str(min_pix_num') repmat(' pix period', [4, 1])]})

%%

f_obj = f_tube/obj_mag;

diameter_bfp1 = 2 * f_obj * obj_NA/obj_RI; % assuming objectives have been optimized to follow paraxial approx

diameter_bfp2 = 2 * f_obj * tan(asin(obj_NA/obj_RI)); % more like reality


NA_effm = SLMm_mm * M_phase / diameter_bfp2; 
NA_effn = SLMn_mm * M_phase / diameter_bfp2; 

%% compute defocus efficiency for microscope
num_mic = numel(microscope_NA_eff);

%resolution_lateral = (0.61*wavelength)./(NA_eff);
%resolution_axial = (2*wavelength)./(NA_eff.^2);

resolution_x_z = zeros(num_mic, 2);
lateral_edge_efficiency = zeros(num_mic, 1);
lateral_eff_lim_dist = zeros(num_mic, 1);
axial_eff_lim_dist = zeros(num_mic, 1);
for n_mic = 1:num_mic
    % lateral
    resolution_x_z(n_mic, 1) = (0.61*wavelength)./(microscope_NA_eff(n_mic));
    % axial
    resolution_x_z(n_mic, 2) = (2*wavelength)./(microscope_NA_eff(n_mic).^2);
    
    % define phase for lateral
    rho = linspace(0, 1, microscope_SLM_pixels(n_mic)/2);
    sin_alpha = microscope_NA_eff(n_mic)/obj_RI;
    phase = obj_RI * k * sin_alpha * rho;
    phase = phase - min(phase);
    
    % compute eff at FOV edge
    x_edge = microscope_FOV_size(n_mic)/2;

    lateral_edge_efficiency(n_mic) = f_comp_lat_eff(phase, x_edge);
    
    % 95% efficiency x loc
    temp_eff = 1;
    temp_x = 0;
    
    while temp_eff > eff_lim
        temp_x = temp_x + d_x;
        temp_eff = f_comp_lat_eff(phase, temp_x);
    end
    lateral_eff_lim_dist(n_mic) = temp_x - d_x;
    
    % 95% eff of defocus
    phase = obj_RI * k * sqrt(1 - rho.^2 * sin_alpha^2);
    phase = phase - min(phase);
    temp_eff = 1;
    temp_z = 0;
    
    while temp_eff>eff_lim
        temp_z = temp_z + d_z;
        temp_eff = f_comp_defocus_eff(phase, temp_z, microscope_SLM_pixels(n_mic));
    end
    axial_eff_lim_dist(n_mic) = temp_z - d_z;
end

%% eff NA vs resolution

NA_eff = 0.1:0.01:1.05;

resolution_lateral = (0.61*wavelength)./(NA_eff);
resolution_axial = (2*wavelength)./(NA_eff.^2);

figure; hold on;
plot(resolution_lateral/1e-6, NA_eff, 'Linewidth', 2);
plot(resolution_axial/1e-6, NA_eff, 'Linewidth', 2);
if add_scope_params
    for n_mic = 1:num_mic
        plot(resolution_x_z(n_mic,1)/1e-6, microscope_NA_eff(n_mic), [microscope_symbols{n_mic} microscope_colors{n_mic}], 'MarkerSize' ,7);
    end
    for n_mic = 1:num_mic
        plot(resolution_x_z(n_mic,2)/1e-6, microscope_NA_eff(n_mic), [microscope_symbols{n_mic} microscope_colors{n_mic}], 'MarkerSize' ,7);
    end
    legend([{'Lateral (x-y)', 'Axial (z)'} microscope_names])
else
    legend({'Lateral (x-y)', 'Axial (z)'});
end
xlabel('Resolution (um)')
ylim([NA_eff(1) NA_eff(end)])
xlim([0, 40])
ylabel('NAeff')

title(sprintf('Rayleigh resolution vs NAeff at %.0fnm', wavelength/1e-9))


%% lateral efficiency

% n*k*x*sin(theta)
% n*k*x*sin(alpha)*rho

num_x = numel(x_range);
num_NA = numel(NA_eff_plot);

rho = linspace(0, 1, SLMm/2);
colors1 = parula(num_NA);

eff_all_x = zeros(num_x ,num_NA);
figure; hold on;
for n_na = 1:num_NA
    sin_alpha = NA_eff_plot(n_na)/obj_RI;
    phase = obj_RI * k * sin_alpha * rho;
    phase = phase - min(phase);

    for n_x = 1:num_x
        eff_all_x(n_x, n_na) = f_comp_lat_eff(phase, x_range(n_x));
    end
    plot(x_range/1e-6, eff_all_x(:,n_na), 'color', colors1(n_na,:), 'Linewidth', 2)
end
if add_scope_params
    for n_mic = 1:num_mic
        plot(microscope_FOV_size(n_mic)/2/1e-6, lateral_edge_efficiency(n_mic), [microscope_symbols{n_mic} microscope_colors{n_mic}], 'MarkerSize' ,7);
    end
    n_mic = 1;
    line([microscope_FOV_size(n_mic)/2/1e-6, microscope_FOV_size(n_mic)/2/1e-6], [0 1], 'color', [0.5 0.5 0.5], 'Linestyle', '--')
    t = text(microscope_FOV_size(n_mic)/2/1e-6+6, 0.01, '25X FOV edge');
    t.Rotation = 90;
    t.FontSize = 8;
    n_mic = 3;
    line([microscope_FOV_size(n_mic)/2/1e-6, microscope_FOV_size(n_mic)/2/1e-6], [0 1], 'color', [0.5 0.5 0.5], 'Linestyle', '--')
    t = text(microscope_FOV_size(n_mic)/2/1e-6+6, 0.01, '20X FOV edge');
    t.Rotation = 90;
    t.FontSize = 8;
end

legend(num2str(NA_eff_plot'), 'location', 'southwest')
xlabel('Lateral diffraction distance (um)')
ylabel('Diffraction efficiency')
xlim('tight')
ylim([0 1]);
title(sprintf('Lateral diffraction efficiency vs NAeff, for %d pixel SLM', SLMm))

%% max lateral efficiency 

num_pat = numel(pat_pixels);
num_NA = numel(NA_eff_comp);

colors1 = parula(num_pat+1);

dist_all = zeros(num_NA, num_pat);
figure; hold on;
for n_pat = 1:num_pat
    rho = linspace(0, 1, pat_pixels(n_pat)/2);
    
    for n_na = 1:num_NA
        sin_alpha = NA_eff_comp(n_na)/obj_RI;
        phase = obj_RI * k * sin_alpha * rho;
        phase = phase - min(phase);

        temp_eff = 1;
        temp_x = 0;
        
        while temp_eff > eff_lim
            temp_x = temp_x + d_x;
            temp_eff =  f_comp_lat_eff(phase, temp_x);
        end
        dist_all(n_na, n_pat) = temp_x - d_x;
    end
    plot(dist_all(:, n_pat)/1e-6, NA_eff_comp, 'Linewidth', 2, 'color', colors1(n_pat,:));
end

if add_scope_params
    for n_mic = 1:num_mic
        plot(lateral_eff_lim_dist(n_mic)/1e-6, microscope_NA_eff(n_mic), [microscope_symbols{n_mic} microscope_colors{n_mic}], 'MarkerSize' ,7);
    end
%     n_mic = 1;
%     line([microscope_FOV_size(n_mic)/2/1e-6, microscope_FOV_size(n_mic)/2/1e-6], [0 1], 'color', [0.5 0.5 0.5], 'Linestyle', '--')
%     t = text(microscope_FOV_size(n_mic)/2/1e-6+6, 0.01, '25X FOV edge');
%     t.Rotation = 90;
%     t.FontSize = 8;
%     n_mic = 3;
%     line([microscope_FOV_size(n_mic)/2/1e-6, microscope_FOV_size(n_mic)/2/1e-6], [0 1], 'color', [0.5 0.5 0.5], 'Linestyle', '--')
%     t = text(microscope_FOV_size(n_mic)/2/1e-6+6, 0.01, '20X FOV edge');
%     t.Rotation = 90;
%     t.FontSize = 8;
end
legend(num2str(pat_pixels'))
xlabel(sprintf('Lateral diffraction distance (um)'))
ylabel('NAeff')
ylim([NA_eff_comp(1) NA_eff_comp(end)]);
xlim([0 max_dist])
title(sprintf('Distance at %d%% lateral diffraction efficiency vs NAeff', eff_lim*100))


%% defocus efficiency for param range

% n*k*z*cos(theta)
% n*k*z*sqrt(1-sin(theta)^2)
% sin(theta) = rho*sin(alpha)

num_z = numel(z_range);
num_NA = numel(NA_eff_plot);

rho = linspace(0, 1, SLMm/2);
colors1 = parula(num_NA);

eff_all = zeros(num_z,num_NA);
figure; hold on;
for n_na = 1:num_NA
    sin_alpha = NA_eff_plot(n_na)/obj_RI;
    phase = obj_RI * k * sqrt(1 - rho.^2 * sin_alpha^2);
    phase = phase - min(phase);

    for n_z = 1:num_z
        eff_all(n_z, n_na) = f_comp_defocus_eff(phase, z_range(n_z), SLMm);
    end
    plot(z_range/1e-6, eff_all(:,n_na), 'color', colors1(n_na,:), 'Linewidth', 2)
end
legend(num2str(NA_eff_plot'), 'location', 'southwest')
xlabel('Defocus distance (um)')
ylabel('Diffraction efficiency')
xlim('tight')
ylim([0 1]);
title(sprintf('Defocus efficiency vs NAeff, for %d pixel SLM', SLMm))

%% max defocus dist

num_pat = numel(pat_pixels);
num_NA = numel(NA_eff_comp);

colors1 = parula(num_pat+1);

dist_all = zeros(num_NA, num_pat);
figure;  hold on
for n_pat = 1:num_pat
    rho = linspace(0, 1, pat_pixels(n_pat)/2);
    
    for n_na = 1:num_NA
        sin_alpha = NA_eff_comp(n_na)/obj_RI;
        phase = obj_RI * k * sqrt(1 - rho.^2 * sin_alpha^2);
        phase = phase - min(phase);

        temp_eff = 1;
        temp_z = 0;
        
        while temp_eff>eff_lim
            temp_z = temp_z + d_z;
            temp_eff = f_comp_defocus_eff(phase, temp_z, pat_pixels(n_pat));
        end
        dist_all(n_na, n_pat) = temp_z - d_z;
    end
    plot(dist_all(:, n_pat)/1e-6, NA_eff_comp, 'Linewidth', 2, 'color', colors1(n_pat,:));
end
if add_scope_params
    for n_mic = 1:num_mic
        plot(axial_eff_lim_dist(n_mic)/1e-6, microscope_NA_eff(n_mic), [microscope_symbols{n_mic} microscope_colors{n_mic}], 'MarkerSize' ,7);
    end
%     n_mic = 1;
%     line([microscope_FOV_size(n_mic)/2/1e-6, microscope_FOV_size(n_mic)/2/1e-6], [0 1], 'color', [0.5 0.5 0.5], 'Linestyle', '--')
%     t = text(microscope_FOV_size(n_mic)/2/1e-6+6, 0.01, '25X FOV edge');
%     t.Rotation = 90;
%     t.FontSize = 8;
%     n_mic = 3;
%     line([microscope_FOV_size(n_mic)/2/1e-6, microscope_FOV_size(n_mic)/2/1e-6], [0 1], 'color', [0.5 0.5 0.5], 'Linestyle', '--')
%     t = text(microscope_FOV_size(n_mic)/2/1e-6+6, 0.01, '20X FOV edge');
%     t.Rotation = 90;
%     t.FontSize = 8;
end
legend(num2str(pat_pixels'))
xlabel(sprintf('Defocus distance (um)'))
ylabel('NAeff')
ylim([NA_eff_comp(1) NA_eff_comp(end)]);
xlim([0 max_dist])
title(sprintf('Distance at %d%% defocus efficiency vs NAeff', eff_lim*100))




%% zernike orthogonality
% 
% d_rho = 0.000001;
% x_loc = 0:d_rho:1;
% rho = abs(x_loc);
% 
% y1 = rho.^2;
% 
% y2 = rho.^4;
% 
% zer02 = sqrt(3) * (2*rho.^2 - 1); % 
% 
% zer04 = sqrt(5) * (6*rho.^4 - 6*rho.^2 + 1); % 
% 
% zer06 = sqrt(15) * (20*rho.^6- 30*rho.^4+ 12*rho.^2- 1);
% 
% figure; hold on;
% plot(x_loc, zer02)
% plot(x_loc, zer04)
% plot(x_loc, zer06)
% 
% 2*pi*sum(zer02.*rho*d_rho)
% 
% 
% figure;
% plot(rho, zer04 - zer02)
% 
% 
% min(zer02)
% 
% max(zer02)


%% defoucs





















