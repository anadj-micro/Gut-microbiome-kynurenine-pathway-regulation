%% making tables of counts for each taxonomical level

% read asv tables
asv = readtable('tblASVcounts_human_filter.csv', 'ReadVariableNames', 1, 'Delimiter', ',');
taxonomy = readtable("tblASVtaxonomy_silva132_v4v5_filter.csv");

% merge two tables
asv1 = outerjoin(asv,taxonomy);

% first for every entry with NA in classification, replace it with the
% classification from the previous level (NA genus, search family and put 'unclassified_FAMILY)

for i = 1:size(asv1,1)
    asv = asv1.ASV_asv{i};
    classificationGenus = asv1.Genus{i};
    if strcmp(classificationGenus, '<not present>') == 1
        classificationFamily = asv1.Family{i};
        if strcmp(classificationFamily, '<not present>') ~= 1
            s =  cellstr(strcat('unclassified', '_', classificationFamily, '_', asv));
            asv1(i, 'Genus') = s;
        else
            classificationOrder = asv1.Order{i};
            if strcmp(classificationOrder, '<not present>') ~= 1
                s =  cellstr(strcat('unclassified', '_', classificationOrder, '_', asv));
                asv1(i, 'Genus') = s;
                asv1(i, 'Family') = s;
            else
                classificationClass = asv1.Class{i};
                if strcmp(classificationClass, '<not present>') ~= 1
                    s =  cellstr(strcat('unclassified', '_', classificationClass, '_', asv));
                    asv1(i, 'Genus') = s;
                    asv1(i, 'Family') = s;
                    asv1(i, 'Order') = s;
                else
                    classificationPhylum = asv1.Phylum{i};
                    if strcmp(classificationPhylum, '<not present>') ~= 1
                        s =  cellstr(strcat('unclassified', '_', classificationPhylum, '_', asv));
                        asv1(i, 'Genus') = s;
                        asv1(i, 'Family') = s;
                        asv1(i, 'Order') = s;
                        asv1(i, 'Class') = s;
                    else
                        asv1(i, 'Genus') = cellstr('unclassified');
                        asv1(i, 'Family') = cellstr('unclassified');
                        asv1(i, 'Order') = cellstr('unclassified');
                        asv1(i, 'Class') = cellstr('unclassified');
                        asv1(i, 'Phylum') = cellstr('unclassified');
                    end
                end
            end
        end
    end
end

% now parse tables for each taxonomic level
taxLevels = asv1.Properties.VariableNames(6:11);
for i = 1:size(taxLevels,2)
    taxLevel = taxLevels{i};
    uniqueMembers = unique(asv1(:,taxLevel));
    subasv1 = asv1(:,{'SampleID', 'Count', taxLevel});
    asv2 = unstack(subasv1, 'Count', taxLevel);
    asv3 = table2array(asv2(:,2:end));
    asv3(isnan(asv3)) = 0;
    asv2{:,2:end} = asv3;
    % get the same but in relative abundances
    sumCounts = sum(asv2{:,2:end},2);
    asv2Rel = asv2;
    asv2Rel{:,2:end} = (asv2{:,2:end}*100)./sumCounts;
    writetable(asv2, strcat(taxLevel, '_dada2.txt'), 'Delimiter', '\t');
    writetable(asv2Rel, strcat(taxLevel, '_percentage_dada2.txt'), 'Delimiter', '\t');
end

% make ASV wide table as well

asv2 = unstack(asv1(:,1:3),'Count', 'ASV_asv');
asv3 = table2array(asv2(:,2:end));
asv3(isnan(asv3)) = 0;
asv2{:,2:end} = asv3;
writetable(asv2, 'ASV_unstacked.txt', 'Delimiter', '\t');

% get relative abundances
asvRel = asv2;
asvRel{:,2:end} = (asv2{:,2:end}*100)./sum(asv2{:,2:end},2);
abundanceSum = sum(asvRel{:, 2:end}, 2)
writetable(asvRel, 'ASV_unstacked_percentage.txt', 'Delimiter', '\t')

% get relative abundances
asvRel = asv2;
asvRel{:,2:end} = asv2{:,2:end}./sum(asv2{:,2:end},2);
abundanceSum = sum(asvRel{:, 2:end}, 2)
writetable(asvRel, 'ASV_unstacked_relabund.txt', 'Delimiter', '\t')


