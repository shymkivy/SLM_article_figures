close all;
clear;
fnames = 'zernike_scan_data_4_8_22_3h_11m.mat';

fpath = 'C:\Users\ys2605\Desktop\stuff\SLM_GUI\SLM_outputs\AO_outputs';

data = load([fpath '\' fnames]);

% 8, 10, 12, 13
mode_data = data.ao_params.mode_data_all{1};

mode = [8, 10];

num_modes = numel(mode);

weights1 = cell(num_modes, 1);
intens1 = cell(num_modes, 1);
intens_mean1 = cell(num_modes, 1);
coef_abs = zeros(num_modes, 3);
for n_mode = 1:num_modes
    
    mode_data2 = mode_data([mode_data.mode] == mode(n_mode));

    [~, idx1] = sort([mode_data2.weight]);

    mode_data2 = mode_data2(idx1);

    weights = mean(reshape([mode_data2.weight], 2, []),1);

    intens = reshape([mode_data2.intensity_sm], 2, []);
    
    weights1{n_mode} = weights;
    intens1{n_mode} = intens;
    intens_mean1{n_mode} = mean(intens,1);
    
    yf = fit(weights', mean(intens,1)', 'gauss1');
    x_fit = min(weights):0.1:max(weights);
    y_fit = yf(x_fit);
    
    coef_abs(n_mode,:) = [yf.a1, yf.b1, yf.c1];
    
    figure; hold on;
    plot(weights', intens'); 
    plot(x_fit, y_fit); 
    title(sprintf('mode %d', mode(n_mode)));
    
end

num_w = numel(weights1{1});

surf_int = ones(num_w, num_w);

surf_int = intens_mean1{1} .* intens_mean1{2}';

inc_fac = 3;
x_fit2 = linspace(min(weights)*inc_fac, max(weights)*inc_fac, 30);
[X, Y] = meshgrid(x_fit2, x_fit2);

surf_int1 = 1*exp(-((X-coef_abs(1,2))/coef_abs(1,3)).^2);  % coef_abs(1,1)

surf_int2 = 1*exp(-((Y-coef_abs(2,2))/coef_abs(2,3)).^2); % coef_abs(2,1)
        
figure; 
h = surf(x_fit2, x_fit2, surf_int1.*surf_int2);
xlabel(sprintf('Z_2^-^2 weight'))
ylabel(sprintf('Z_2^2 weight'))
set(h,'edgecolor','k');
set(h,'EdgeAlpha',0.5);
%set(h,'facecolor','none');
%set(h,'FaceAlpha',0.5);
colormap(parula)
grid off  
shading interp
