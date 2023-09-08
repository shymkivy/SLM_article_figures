
clear;
close all


save_path = 'C:\Users\ys2605\Desktop\stuff\papers\SLM_microscope\fig_data\';

%fname = '8_25_23_disk_est';
fname = '8_24_23_point_est';
%fname = '8_25_23_disk_est2';
load([save_path, '\', fname, '.mat']);



colors1 = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0.9290 0.6940 0.1250], [0.3010 0.7450 0.9330]};

colors1 = {[0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]};

colors1 = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]};

colors1 = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0.9290 0.6940 0.1250], [0.3010 0.7450 0.9330]};

%%

plot_idx = logical([1, 1, 0, 1, 1, 1, 0]);

plot_idx = logical([1, 1, 1, 1, 1, 0]);

plot_idx = logical([1, 1, 1, 1]);
alpha1 = 0.1;


%%

algorithms = data_all.params.algorithms;
legend1 = data_all.params.legend1;
num_pts_sample = data_all.params.num_pts_sample;
num_reps = data_all.params.num_reps;
intensity_range = data_all.params.intensity_range;

mag_err_all = data_all.mag_err_all;
peak_eff_all = data_all.peak_eff_all;
zero_ord_all = data_all.zero_ord_all;
compute_dur_all = data_all.compute_dur_all;
peak_abs_err_all = data_all.peak_abs_err_all;


%%

num_alg = numel(algorithms);
num_intens = numel(intensity_range);


for n_leg = 1:num_alg
    if strcmpi(legend1{n_leg}(end-2:end), '_LW')
        legend1{n_leg} = legend1{n_leg}(1:end-3);
    end
end

legend1 = {'PS * global_GS disk',...
           'OPS * global_GS disk',...
           'OP + NOVO_CGH_VarIEuclid disk',...
           'random superposition',...
           'global-GS',...
           'NOVO-CGH VarI',...
           'NOVO-CGH VarIEuclid'};%, 'global_GS_LW'};


legend1 = {'PS',...
           'OPS',...
           'random superposition',...
           'global-GS',...
           'NOVO-CGH VarI',...
           'NOVO-CGH VarIEuclid'};%, 'global_GS_LW'};

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
title(sprintf('Mean absolute error of point magnitudes; %d reps; %d intens\n', num_reps, num_intens));
xlim([0, num_pts_sample(end)])
ylim([0, 5])
%f2.Children(1).FontSize = 10;


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
title(sprintf('Peak efficiency; %d reps; %d intens\n', num_reps, num_intens));
%f3.Children(1).FontSize = 10;
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
title(sprintf('Intensity sent to zero order; %d reps; %d intens\n', num_reps, num_intens));
%f4.Children(1).FontSize = 10;



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
title(sprintf('Computation duration; %d reps; %d intens\n', num_reps, num_intens));
%f5.Children(1).FontSize = 10;

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
title(sprintf('Intensity error within disk; %d reps; %d intens\n', num_reps, num_intens));
%f6.Children(1).FontSize = 10;
%f1.Children(2).FontSize = 10
%f1.Children(2).XScale = 'linear'