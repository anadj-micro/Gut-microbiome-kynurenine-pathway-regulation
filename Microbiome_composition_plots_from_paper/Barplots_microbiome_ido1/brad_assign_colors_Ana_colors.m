function [taxonomy, color_spectrum_dict] = brad_assign_colors_Ana_colors(taxonomy, relative_abundance_matrix)

assigned_colors = cell(size(taxonomy,1),1); % colors for each asv
brad_order_colors = []; % colors based on brad_color
is_tax_assigned_a_color = false(size(assigned_colors)); % which asv has yet to be assigned?

total_asv_abundance_in_all_samples = sum(relative_abundance_matrix,1);
max_asv_abundance_in_all_samples = max(relative_abundance_matrix,[],1);

%% special color settings
fixed_color_taxa = {'Enterococcus', 'Streptococcus', 'Staphylococcus', 'Erysipelotrichales', 'Bacilli', 'Clostridia', ...
    'Bacteroidota', 'Proteobacteria', 'Verrucomicrobiota', 'Others'};

fixed_color_tax_category = {'Genus', 'Genus', 'Genus', 'Order', 'Class', 'Class',...
    'Phylum', 'Phylum', 'Phylum', ''};

% Green, Yellow-green, Yellow, Orange, Blue, Brown, Teal, Red, Purple, Grey
color_lower_limit = {'#006026', '#b9db3f', '#fcf628', '#ffa735', '#708df4', '#d6bdad','#1efff3', '#ff3838', '#dd05ff', '#e2e2e2'};
color_upper_limit = {'#3ce83f', '#8ca530', '#d6d120', '#ed9017', '#012191', '#685c54','#16ddd3', '#af0000', '#862296', '#8e8e8e'};

% % Green, Yellow-green, Yellow, Blue, Brown, Orange, Teal, Red, Purple, Grey
% color_lower_limit = {'#006026', '#b9db3f', '#fcf628', '#708df4', '#d6bdad', '#ffa735', '#1efff3', '#ff3838', '#dd05ff', '#e2e2e2'};
% color_upper_limit = {'#3ce83f', '#8ca530', '#d6d120', '#012191', '#685c54', '#ed9017', '#16ddd3', '#af0000', '#862296', '#8e8e8e'};

color_spectrum_dict = struct(); % map from fixed-color taxonomy to color spectrum generated for all ASVs in this taxonomy
for i=1:length(fixed_color_taxa)
    curr_taxa = fixed_color_taxa{i};
    switch fixed_color_tax_category{i}
        case 'Genus'
            asv_to_assign_a_color = strcmpi(taxonomy.genus,curr_taxa); % ignore case
        case 'Family'
            asv_to_assign_a_color = strcmpi(taxonomy.family,curr_taxa); % ignore case
        case 'Order'
            asv_to_assign_a_color = strcmpi(taxonomy.order,curr_taxa); % ignore case
        case 'Class'
            asv_to_assign_a_color = strcmpi(taxonomy.class,curr_taxa); % ignore case
        case 'Phylum'
            asv_to_assign_a_color = strcmpi(taxonomy.phylum,curr_taxa); % ignore case
        case 'Kingdom'
            asv_to_assign_a_color = strcmpi(taxonomy.kingdom,curr_taxa); % ignore case
        case ''
            % do nothing
        otherwise
            error('%s not recognized\n', fixed_color_tax_category{i});
    end
    
    % ASVs must have not been assigned a color before
    if (fixed_color_tax_category{i})
        asv_to_assign_a_color = asv_to_assign_a_color & ~is_tax_assigned_a_color;
        is_tax_assigned_a_color(asv_to_assign_a_color) = true;
    else
        asv_to_assign_a_color = ~is_tax_assigned_a_color;
    end
    
    % generate a spectrum of colors between color_lower_limit and
    % color_upper_limit
    % nnz gives number of nonzero matrix elements.
    if (nnz(asv_to_assign_a_color)>0)
        [assigned_colors, color_spectrum] = brad_generate_colors_taxonomic_group(...
            assigned_colors,...
            asv_to_assign_a_color,...
            color_lower_limit{i},...
            color_upper_limit{i},...
            total_asv_abundance_in_all_samples,...
            max_asv_abundance_in_all_samples);
        color_spectrum_dict.(fixed_color_taxa{i}) = color_spectrum;
        brad_order_colors = [brad_order_colors; color_spectrum];
    else
        color_spectrum_dict.(fixed_color_taxa{i}) = {''};
    end
end

taxonomy.HexColor = assigned_colors;
taxonomy.ColorOrder = zeros(height(taxonomy),1);
for k = 1:length(brad_order_colors)
    taxonomy.ColorOrder(find(strcmp(taxonomy.HexColor, brad_order_colors{k}))) = k;
end

end