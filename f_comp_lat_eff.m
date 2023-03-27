function eff = f_comp_lat_eff(phase, x)
%pix_levels = floor((abs(z_range(n_z)) * phase)/(2*pi));
phase1 = x * phase;
phase1 = phase1 - min(phase1);
pix_levels = floor(floor((phase1)/(pi))/2);

[~,level_changes, ~] = unique(pix_levels, 'last');
level_counts = diff([0; level_changes]);

level_eff = (sin(pi./level_counts)./(pi./level_counts)).^2;
level_frac_areas = level_counts/sum(level_counts);

eff = sum(level_eff.*level_frac_areas);

end