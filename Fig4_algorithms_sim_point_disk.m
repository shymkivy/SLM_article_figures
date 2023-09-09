close all
clear

%%
% need to download the slm gui and provide path
slm_dir = 'C:\Users\ys2605\Desktop\stuff\SLM_GUI\SLM_GUI';
addpath(genpath([slm_dir '\SLM_GUI_funcions']))

% make a save path first
save_path = 'C:\Users\ys2605\Desktop\stuff\papers\SLM_microscope\fig_data\';
save_name = '8_25_23_disk_est3';

view_z = 0;
point_rad_um = 5;        % area to integrate intensity over, make same as disk for disks
I_est_2P = false;       % check if you want to estimate squared intensity in simulation
plot_stuff = 0;
make_disk = 1;          % set to 1 if want to simulate disks
disk_rad_um = 5;        % radius of disks to simulate


num_pts_sample = [2, 4, 10, 20, 50, 100, 200]; % points to simulate %
num_reps = 30; % repeats per group

x_samp = -250:250;      % coordinates to sample from
y_samp = -250:250;
z_samp = -250:25:250;

intensity_range = linspace(0.1,1,100);  % intensities to sample from for each point
%intensity_range = 0;                      % for same intensity everywhere

ignore_coords = -5:5;  % ignore x-y region because of zero order measurement


% indicate which algorigms to run, their legends

% for disk generation (point holo alg, legend, disk holo alg)
algo_leg_disk = {'superposition',            'phase superpos + globalGS disk',                'global_GS_LW';...
                 'superposition_optim',      'optimized phase superpos + globalGS disk',      'global_GS_LW';...
                 'superposition',            'phase superpos + NOVO_CGH_VarIEuclid disk',     'NOVO_CGH_VarIEuclid_LW';...      
                 'superposition_LW',         'random superposition',                          '';...
                 'global_GS_LW',             'global_GS_LW',                                  '';...
                 'NOVO_CGH_VarI_LW',         'NOVO_CGH_VarI_LW',                              '';...
                 'NOVO_CGH_VarIEuclid_LW',   'NOVO_CGH_VarIEuclid_LW',                        '';...
                 };

% for point generation (point holo alg, legend)
algo_leg =      {'superposition',             'phase superpos';...
                 'superposition_optim',       'optimized phase superpos';...
                 'superposition_LW',          'random superposition';...      
                 'global_GS_LW',              'global GS';...
                 'NOVO_CGH_VarI_LW',          'NOVO_CGH_VarI';...
                 'NOVO_CGH_VarIEuclid_LW',    'NOVO_CGH_VarIEuclid';...
                 };

% index of which algorithm inputs to plot
ignore_in_plot_idx = [];
%ignore_in_plot_idx = [2, 4];

alpha1 = 0.1;       % for transparancy level

colors1 = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0.9290 0.6940 0.1250], [0.3010 0.7450 0.9330]};

colors1 = {[0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]};


%% SLM parameters

reg1.SLMm = 1152;
reg1.SLMn = 960;
reg1.wavelength = 940;
reg1.effective_NA = 0.6320;
reg1.objective_RI = 1.33;
reg1.phase_diameter = 1152;
reg1.xyz_affine_tf_mat = eye(3);
reg1.xyz_offset = [0, 0, 0];
reg1.zero_outside_phase_diameter = 1;
reg1.holo_mask = true(reg1.SLMm, reg1.SLMn);

if reg1.zero_outside_phase_diameter
    xlm = linspace(-reg1.SLMm/reg1.phase_diameter, reg1.SLMm/reg1.phase_diameter, reg1.SLMm);
    xln = linspace(-reg1.SLMn/reg1.phase_diameter, reg1.SLMn/reg1.phase_diameter, reg1.SLMn);
    [fX, fY] = meshgrid(xln, xlm);
    [~, RHO] = cart2pol(fX, fY);
    reg1.holo_mask(RHO>1) = 0;
end

coord.xyzp = [0, 0, 0];
coord.I_targ = 1;
coord.I_targ1P = 1;
coord.W_est = 1;

% make compatible with GUI functions
app.MakedisksCheckBox.Value = make_disk;
app.DiskradiusumEditField.Value = disk_rad_um;
app.UsegaussianbeamampCheckBox.Value = 1;
app.pointsizeumEditField.Value = point_rad_um;
app.I_estI22PCheckBox.Value = I_est_2P;
apply_ao = 0;


%%

if make_disk
    algorithms = algo_leg_disk(:,1);
    legend1 = algo_leg_disk(:,2);
    phase_synthesis_disk_method = algo_leg_disk(:,3);
    app.pointsizeumEditField.Value = disk_rad_um;
else
    algorithms = algo_leg(:,1);
    legend1 = algo_leg(:,2);
    phase_synthesis_disk_method = cell(size(algo_leg,1));
end
plot_idx = true(numel(algorithms),1);

plot_idx(ignore_in_plot_idx) = 0;

params.reg1 = reg1;
params.algorithms = algorithms;
params.legend1 = legend1;
params.num_pts_sample = num_pts_sample;
params.num_reps = num_reps;
params.x_samp = x_samp;
params.y_samp = y_samp;
params.z_all = z_samp;
params.intensity_range = intensity_range;
params.ignore_coords = ignore_coords;

%%

num_samp = numel(num_pts_sample);
num_alg = numel(algorithms);
num_intens = numel(intensity_range);
num_z = numel(z_samp);

x_samp(logical(sum(x_samp == ignore_coords',1))) = [];
y_samp(logical(sum(y_samp == ignore_coords',1))) = [];

mag_err_all = cell(num_alg,1);
zero_ord_all = cell(num_alg,1);
compute_dur_all = cell(num_alg,1);
peak_eff_all = cell(num_alg,1);
peak_abs_err_all = cell(num_alg,1);


for n_alg = 1:num_alg
    mag_error = zeros(num_samp, num_reps);
    zero_ord_mag = zeros(num_samp, num_reps);
    compute_duration = zeros(num_samp, num_reps);
    peak_eff = zeros(num_samp, num_reps);
    peak_abs_err = zeros(num_samp, num_reps);
    fprintf('Computing: %s; samp x/%d..', algorithms{n_alg}, num_samp);
    for n_samp = 1:num_samp
        num_pts1 = num_pts_sample(n_samp);
        fprintf(' %d', n_samp);
        for n_rep = 1:num_reps
    
            %x_pts = round(rand(num_pts1,1)*x_span + x_range(1));
            %y_pts = round(rand(num_pts1,1)*y_span + y_range(1));
            x_pts = randsample(x_samp, num_pts1, true)';
            y_pts = randsample(y_samp, num_pts1, true)';
            z_pts = z_samp(randsample(num_z, num_pts1, true))';
            
            intens = intensity_range(randsample(num_intens, num_pts1, true))';

            coord.xyzp = [x_pts, y_pts, z_pts];
            coord.I_targ = intens;
            coord.I_targ1P = intens;
            coord.W_est = sqrt(intens);
            
            tic();
            if strcmpi(algorithms{n_alg}, 'superposition_optim')
                [SLM_phase, holo_phase, SLM_phase_corr, holo_phase_corr, AO_phase] = f_sg_xyz_gen_SLM_phase(app, coord, reg1, apply_ao, 'superposition', phase_synthesis_disk_method{n_alg});
                w_out = f_sg_optimize_phase_w(app, reg1, holo_phase, coord, coord.I_targ, 0);
                coord2 = coord;
                coord2.W_est = w_out.w_final;
                [SLM_phase, holo_phase, SLM_phase_corr, holo_phase_corr, AO_phase] = f_sg_xyz_gen_SLM_phase(app, coord2, reg1, apply_ao, 'superposition', phase_synthesis_disk_method{n_alg});
            else
                [SLM_phase, holo_phase, SLM_phase_corr, holo_phase_corr, AO_phase] = f_sg_xyz_gen_SLM_phase(app, coord, reg1, apply_ao, algorithms{n_alg});
            end
            compute_duration(n_samp, n_rep) = toc();
    
            if reg1.zero_outside_phase_diameter
                SLM_phase(~reg1.holo_mask) = 0;
                SLM_phase_corr(~reg1.holo_mask) = 0;
            end
        
            data_w = f_sg_simulate_intensity(reg1, SLM_phase, coord, point_rad_um, app.UsegaussianbeamampCheckBox.Value, I_est_2P, 0);
            
            peak_abs_err(n_samp, n_rep) = mean(data_w.pt_abs_err);

            mags_corr = data_w.pt_mags./coord.I_targ;

            mean1 = mean(mags_corr);
            mag_error(n_samp, n_rep) = mean(abs(mags_corr - mean1)/mean1*100);
    
            zero_ord_mag(n_samp, n_rep) = data_w.zero_ord_mag/mean(data_w.im_sum)*100;

            peak_eff(n_samp, n_rep) = sum(data_w.pt_mags)/mean(data_w.im_sum)*100;
    
        end
    end
    
    fprintf('\n');
    
    mag_err_all{n_alg} = mag_error;
    zero_ord_all{n_alg} = zero_ord_mag;
    compute_dur_all{n_alg} = compute_duration;
    peak_eff_all{n_alg} = peak_eff;
    peak_abs_err_all{n_alg} = peak_abs_err;
end  
fprintf('Done\n');

data_all.mag_err_all = mag_err_all;
data_all.zero_ord_all = zero_ord_all;
data_all.compute_dur_all = compute_dur_all;
data_all.peak_eff_all = peak_eff_all;
data_all.peak_abs_err_all = peak_abs_err_all;

params.reg1 = reg1;
params.num_pts_sample = num_pts_sample;
params.algorithms = algorithms;
params.legend1 = legend1;
params.num_reps = num_reps;
params.z_all = z_samp;
params.intensity_range = intensity_range;
params.ignore_coords = ignore_coords;
params.x_samp = x_samp;
params.y_samp = y_samp;

data_all.params = params;

%%
f1 = figure; hold on;
pl_all = cell(num_alg,1);
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        mag_error = mag_err_all{n_alg};
        plot(num_pts_sample, 100-mag_error, 'color', [colors1{n_alg}, alpha1])
    end
end
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        mag_error = mag_err_all{n_alg};
        pl_all{n_alg} = plot(num_pts_sample, 100-mean(mag_error,2), 'color', colors1{n_alg}, LineWidth=2);
        errorbar(num_pts_sample, 100-mean(mag_error,2), std(mag_error, [], 2)/sqrt(num_reps-1), 'color', colors1{n_alg}, LineWidth=2)
    end
end
legend([pl_all{plot_idx}], legend1(plot_idx), 'interpreter', 'none', Location='Southwest')
ylabel('100 - Mean abs error (%)');
xlabel('number points simulated');
title(sprintf('\nMean absolute error of point magnitudes; %d reps; %d intens', num_reps, num_intens));
f1.Children(1).FontSize = 10;
ylim([90, 100])
%f1.Children(2).FontSize = 10
%f1.Children(2).XScale = 'linear'


f2 = figure; hold on;
pl_all = cell(num_alg,1);
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        mag_error = mag_err_all{n_alg};
        plot(num_pts_sample, mag_error, 'color', [colors1{n_alg}, alpha1])
    end
end
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        mag_error = mag_err_all{n_alg};
        pl_all{n_alg} = plot(num_pts_sample, mean(mag_error,2), 'color', colors1{n_alg}, LineWidth=2);
        errorbar(num_pts_sample, mean(mag_error,2), std(mag_error, [], 2)/sqrt(num_reps-1), 'color', colors1{n_alg}, LineWidth=2)
    end
end
legend([pl_all{plot_idx}], legend1(plot_idx), 'interpreter', 'none')
ylabel('Mean abs error (%)');
xlabel('number points simulated');
title(sprintf('Mean absolute error of point magnitudes; %d reps; %d intens', num_reps, num_intens));
%ylim([0, 5])
xlim([0, num_pts_sample(end)])
f2.Children(1).FontSize = 10;

f3 = figure; hold on;
pl_all = cell(num_alg,1);
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        peak_eff = peak_eff_all{n_alg};
        plot(num_pts_sample, peak_eff, 'color', [colors1{n_alg}, alpha1])
    end
end
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        peak_eff = peak_eff_all{n_alg};
        pl_all{n_alg} = plot(num_pts_sample, mean(peak_eff,2), 'color', colors1{n_alg}, LineWidth=2);
        errorbar(num_pts_sample, mean(peak_eff,2), std(peak_eff, [], 2)/sqrt(num_reps-1), 'color', colors1{n_alg}, LineWidth=2)
    end
end
legend([pl_all{plot_idx}], legend1(plot_idx), 'interpreter', 'none')
ylabel('Intensity in peaks (%)');
xlabel('number points simulated');
title(sprintf('Peak efficiency; %d reps; %d intens', num_reps, num_intens));
f3.Children(1).FontSize = 10;
%ylim([75, 90])

f4 = figure; hold on;
pl_all = cell(num_alg,1);
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        zero_ord_mag = zero_ord_all{n_alg};
        plot(num_pts_sample, zero_ord_mag, 'color', [colors1{n_alg}, 0.2])
    end
end
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        zero_ord_mag = zero_ord_all{n_alg};
        pl_all{n_alg} = plot(num_pts_sample, mean(zero_ord_mag,2), 'color', colors1{n_alg}, LineWidth=2);
        errorbar(num_pts_sample, mean(zero_ord_mag,2), std(zero_ord_mag, [], 2)/sqrt(num_reps-1), 'color', colors1{n_alg}, LineWidth=2)
    end
end
legend([pl_all{plot_idx}], legend1(plot_idx), 'interpreter', 'none')
ylabel('Percent intensity');
xlabel('number points simulated');
title(sprintf('Intensity sent to zero order; %d reps; %d intens', num_reps, num_intens));
f4.Children(1).FontSize = 10;



f5 = figure; hold on;
pl_all = cell(num_alg,1);
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        compute_dur = compute_dur_all{n_alg};
        plot(num_pts_sample, compute_dur, 'color', [colors1{n_alg}, 0.2])
    end
end
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        compute_dur = compute_dur_all{n_alg};
        pl_all{n_alg} = plot(num_pts_sample, mean(compute_dur,2), 'color', colors1{n_alg}, LineWidth=2);
        errorbar(num_pts_sample, mean(compute_dur,2), std(compute_dur, [], 2)/sqrt(num_reps-1), 'color', colors1{n_alg}, LineWidth=2)
    end
end
legend([pl_all{plot_idx}], legend1(plot_idx), 'interpreter', 'none', Location='Southeast')
ylabel('Tiime (sec)');
xlabel('number points simulated');
title(sprintf('Computation duration; %d reps; %d intens', num_reps, num_intens));
f5.Children(1).FontSize = 10;


f6 = figure; hold on;
pl_all = cell(num_alg,1);
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        peak_abs_err = peak_abs_err_all{n_alg};
        plot(num_pts_sample, peak_abs_err, 'color', [colors1{n_alg}, alpha1])
    end
end
for n_alg = 1:num_alg
    if plot_idx(n_alg)
        peak_abs_err = peak_abs_err_all{n_alg};
        pl_all{n_alg} = plot(num_pts_sample, mean(peak_abs_err,2), 'color', colors1{n_alg}, LineWidth=2);
        errorbar(num_pts_sample, mean(peak_abs_err,2), std(peak_abs_err, [], 2)/sqrt(num_reps-1), 'color', colors1{n_alg}, LineWidth=2)
    end
end
legend([pl_all{plot_idx}], legend1(plot_idx), 'interpreter', 'none', Location='Southwest')
ylabel('Abs error (%)');
xlabel('number points simulated');
title(sprintf('Intensity error within disk; %d reps; %d intens', num_reps, num_intens));
f6.Children(1).FontSize = 10;
%f1.Children(2).FontSize = 10
%f1.Children(2).XScale = 'linear'

%%

save([save_path, '\', save_name], 'data_all', "-v7.3");

