function [assigned_colors,color_spectrum] = brad_generate_colors_taxonomic_group(...
            assigned_colors, asv_to_assign_a_color, color_lower_limit, color_upper_limit,...
            total_abundance_in_all_samples_for_all_asvs, max_asv_abundance_in_all_samples)
       
% maximum number of colors between its lower limit and upper limit (best if odd)
% use the first max_num_colors-1 colors for significant asvs and use
% the last color for non-significant asvs
max_num_colors = 21;  

% assign colors only for significant ASVs (total relative abundance > 5% in
% all samples)
% if the number of significant ASVs is larger than max_num_colors, only
% assign colors to max_num_color ASVs
asv_to_assign_a_color_indices = find(asv_to_assign_a_color);
all_significant_asvs = (max_asv_abundance_in_all_samples > 0.05)';
total_abundance_in_all_samples_for_asvs_to_assign_a_color = total_abundance_in_all_samples_for_all_asvs(asv_to_assign_a_color);
significant_asvs_to_assign_a_color_indices = find(asv_to_assign_a_color & all_significant_asvs);

if length(significant_asvs_to_assign_a_color_indices)>(max_num_colors-1)
    % only color "max_num_colors" significant asvs
    [~,sorted_significant_asv_indices] = sort(total_abundance_in_all_samples_for_all_asvs(significant_asvs_to_assign_a_color_indices), 'descend');
    ordered_significant_asv_indices = significant_asvs_to_assign_a_color_indices(sorted_significant_asv_indices(1:max_num_colors-1));   %ordered color
    leftover_asv_indices = setdiff(asv_to_assign_a_color_indices, ordered_significant_asv_indices);
else
    % color the most abundant "max_num_colors" asvs
    [~,sorted_significant_asv_indices] = sort(total_abundance_in_all_samples_for_all_asvs(asv_to_assign_a_color_indices), 'descend');
    ordered_significant_asv_indices = asv_to_assign_a_color_indices(sorted_significant_asv_indices(1:min(max_num_colors-1,length(total_abundance_in_all_samples_for_asvs_to_assign_a_color))));
    leftover_asv_indices = setdiff(asv_to_assign_a_color_indices, ordered_significant_asv_indices);
end

% generate "max_num_colors-1" colors for significant asvs and one color for the leftover
[significant_colors, non_significant_color] = generate_distinguishablecolors(...
    color_lower_limit,color_upper_limit,length(ordered_significant_asv_indices));
color_spectrum = [significant_colors; non_significant_color];
assigned_colors(ordered_significant_asv_indices) = significant_colors;
assigned_colors(leftover_asv_indices) = non_significant_color;

end

