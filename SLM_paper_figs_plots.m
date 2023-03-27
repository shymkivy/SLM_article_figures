close all
clear


%%
SLMm = 1152;
SLMn = 960; % half width of 1920


SLMm_mm = 10.7;
SLMn_mm = 8.8; % half width 


obj_mag = 25;
obj_NA = 1.05;
% 
% obj_mag = 20;
% obj_NA = .95;

wavelength = 940e-9;

k = 2*pi/wavelength;

obj_RI = 1.33;
f_tube = 180;  % in mm

f_obj = f_tube/obj_mag;



%% diffraction efficiency of blazed gratings vs number of levels

N = 1:100;% number of levels in blazed

nu = (sin(pi./N)./(pi./N)).^2;

figure;
subplot(1,2,1);
plot(N, nu, 'k', 'linewidth', 2); 
axis tight;
title('Blazed grating diffraction efficiency')
ylabel('diffraction efficiency');
xlabel('Number of grating levels per period');

subplot(1,2,2);
semilogx(N, nu, 'k', 'linewidth', 2); 
axis tight;
title('Blazed grating diffraction efficiency')
ylabel('diffraction efficiency');
xlabel('Number of grating levels per period');

figure;
semilogx(N, nu, 'k', 'linewidth', 2); 
axis tight;
title('Blazed grating diffraction efficiency')
ylabel('diffraction efficiency');
xlabel('Number of grating levels per period');


%%

M_phase = 56/55; % magnification of SLM phase to back of objective
M_image = 33/1400; % magnification of points image pattern from right after SLM

wavelength = 9.4e-7;
%SLM_pix_pitch = 9.2e-6;

SLM_pix_pitch = (5:0.2:12)*1e-6;

min_pix_num = [2, 4, 8, 16];
SLM_def_dist_all = zeros(numel(SLM_pix_pitch), numel(min_pix_num));

for n_px = 1:numel(min_pix_num)

    % period*sin(theta) = m*wavelength
    SLM_def_angle = asin(wavelength./(SLM_pix_pitch*min_pix_num(n_px)));

    f_lens = .3;

    SLM_def_dist = f_lens*tan(SLM_def_angle) * M_image;
    
    SLM_def_dist_all(:,n_px) = SLM_def_dist;
end
figure; 
plot(SLM_pix_pitch/1e-6, SLM_def_dist_all./1e-6, 'linewidth', 2)
xlabel('pixel pitch (um)')
ylabel('deflection dist (um)');
title('Blazed period va deflection distance')
legend({[num2str(min_pix_num') repmat(' pix period', [4, 1])]})

%%

diameter_bfp1 = 2 * f_obj * obj_NA/obj_RI; % assuming objectives have been optimized to follow paraxial approx

diameter_bfp2 = 2 * f_obj * tan(asin(obj_NA/obj_RI)); % more like reality


NA_effm = SLMm_mm * M_phase / diameter_bfp2; 
NA_effn = SLMn_mm * M_phase / diameter_bfp2; 
%%

% n*k*z*cos(theta)
% n*k*z*sqrt(1-sin(theta)^2)
% sin(theta) = rho*sin(alpha)

eff_NA = .6;

z_range = (-500:10:500)*1e-6;

num_z = numel(z_range);

rho = linspace(0, 1, SLMm/2);
sin_alpha = eff_NA/obj_RI;

phase = obj_RI * k * sqrt(1 - rho.^2 * sin_alpha^2);
phase = phase - min(phase);

eff_all = zeros(num_z,1);
for n_z = 1:num_z
    %pix_levels = floor((abs(z_range(n_z)) * phase)/(2*pi));
    pix_levels = floor(floor((abs(z_range(n_z)) * phase)/(pi))/2);
    

    [C,level_changes, ic] = unique(pix_levels, 'last');
    level_counts = -diff([level_changes; 0]);
    
    level_eff = (sin(pi./level_counts)./(pi./level_counts)).^2;
    
    level_rad = level_changes/(SLMm/2);
    frac_areas = pi*level_rad.^2/pi;
    level_frac_areas = -diff([frac_areas; 0]);
    
    eff_all(n_z) = sum(level_eff.*level_frac_areas);
end

figure; 
plot(z_range, eff_all)

figure; 
plot(angle(exp(1i*(phase-pi))))



%% zernike orthogonality

d_rho = 0.000001;
x_loc = 0:d_rho:1;
rho = abs(x_loc);

y1 = rho.^2;

y2 = rho.^4;

zer02 = sqrt(3) * (2*rho.^2 - 1); % 

zer04 = sqrt(5) * (6*rho.^4 - 6*rho.^2 + 1); % 

zer06 = sqrt(15) * (20*rho.^6- 30*rho.^4+ 12*rho.^2- 1);

figure; hold on;
plot(x_loc, zer02)
plot(x_loc, zer04)
plot(x_loc, zer06)

2*pi*sum(zer02.*rho*d_rho)


figure;
plot(rho, zer04 - zer02)


min(zer02)

max(zer02)


%% defoucs





















