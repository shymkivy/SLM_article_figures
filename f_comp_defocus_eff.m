function eff = f_comp_defocus_eff(phase, z, SLMm)

%pix_levels = floor((abs(z_range(n_z)) * phase)/(2*pi));
pix_levels = floor(floor((abs(z) * phase)/(pi))/2);


[~,level_changes, ~] = unique(pix_levels, 'last');
level_counts = -diff([level_changes; 0]);

level_eff = (sin(pi./level_counts)./(pi./level_counts)).^2;

level_rad = level_changes/(SLMm/2);
frac_areas = pi*level_rad.^2/pi;
level_frac_areas = -diff([frac_areas; 0]);

eff = sum(level_eff.*level_frac_areas);

end