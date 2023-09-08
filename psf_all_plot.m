clear;
close all;


fpath = 'C:\Users\ys2605\Desktop\stuff\data\PSF_data\SLM_25x_AO\combined\';

fnames = {
          'psf_data_Original path 5_20_23_2023_6_8.mat';...
          'psf_data_ETL-10-30-4f (fast) 5_18_23_2023_6_9.mat';...
          'psf_data_ETL-16-40-Obj (slow) 5_19_23_2023_6_9.mat';...
          'psf_data_SLM no AO 4_10_23_2023_5_21.mat';...
          'psf_data_SLM AO 5_21_23_2023_6_9.mat';...
          };
      
      
num_files = numel(fnames);

z_loc_all = cell(num_files,1);
fwhm_all = cell(num_files,1);
name_tag_all = cell(num_files,1);
for n_f = 1:num_files
    data = load([fpath '\' fnames{n_f}]);
    z_loc_all{n_f} = data.data_fwhm.z_loc_unique;
    fwhm_all{n_f} = data.data_fwhm.fwhm;
    name_tag_all{n_f} = data.data_fwhm.description;
end


colors1 = parula(num_files);

colors1 = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;...
           0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];


d_tags = {'y', 'x', 'z'};
for n_d = 1:3
    figure; hold on;
    pl = cell(num_files,1);
    for n_f = num_files:-1:1
        num_pts = size(fwhm_all{n_f},1);
        means = zeros(num_pts, 3);
        sems = zeros(num_pts, 3); 
        for n_pt = 1:num_pts
            means(n_pt,:) = mean(fwhm_all{n_f}{n_pt},1);
            sems(n_pt,:) = std(fwhm_all{n_f}{n_pt}, [], 1)/max(sqrt(num_pts-1),1);
        end
        pl{n_f} = plot(z_loc_all{n_f}, means(:,n_d), '.-', 'markersize', 20, 'color', colors1(n_f,:), 'linewidth', 2);
        for n_pt = 1:num_pts
            plot(z_loc_all{n_f}(n_pt), fwhm_all{n_f}{n_pt}(:,n_d), 'o', 'color', colors1(n_f,:));
        end
        errorbar(z_loc_all{n_f}, means(:,n_d), sems(:,n_d), 'color', colors1(n_f,:), 'linewidth', 1);
    end
    legend([pl{:}], name_tag_all')
    title(sprintf('%s fwhm', d_tags{n_d}));
    xlabel('z defocus (um)');
    ylabel('FWHM (um)');
end
