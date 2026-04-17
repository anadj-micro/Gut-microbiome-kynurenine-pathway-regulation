%% making tables of counts for each taxonomical level

% read asv tables
asv = readtable('tblASVcounts.csv', 'ReadVariableNames', 1, 'Delimiter', ',');
taxonomy = readtable('tblASVtaxonomy_dada2.csv');

% merge two tables
asv1 = outerjoin(asv,taxonomy);

% first for every entry with NA in classification, replace it with the
% classification from the previous level (NA genus, search family and put 'unclassified_FAMILY)

for i = 1:size(asv1,1)
    asv = asv1.ASV_asv{i};
    classificationGenus = asv1.genus{i};
    if strcmp(classificationGenus, '<not present>') == 1
        classificationFamily = asv1.family{i};
        if strcmp(classificationFamily, '<not present>') ~= 1
            s =  cellstr(strcat('unclassified', '_', classificationFamily, '_', asv));
            asv1(i, 'genus') = s;
        else
            classificationOrder = asv1.order{i};
            if strcmp(classificationOrder, '<not present>') ~= 1
                s =  cellstr(strcat('unclassified', '_', classificationOrder, '_', asv));
                asv1(i, 'genus') = s;
                asv1(i, 'family') = s;
            else
                classificationClass = asv1.class{i};
                if strcmp(classificationClass, '<not present>') ~= 1
                    s =  cellstr(strcat('unclassified', '_', classificationClass, '_', asv));
                    asv1(i, 'genus') = s;
                    asv1(i, 'family') = s;
                    asv1(i, 'order') = s;
                else
                    classificationPhylum = asv1.phylum{i};
                    if strcmp(classificationPhylum, '<not present>') ~= 1
                        s =  cellstr(strcat('unclassified', '_', classificationPhylum, '_', asv));
                        asv1(i, 'genus') = s;
                        asv1(i, 'family') = s;
                        asv1(i, 'order') = s;
                        asv1(i, 'class') = s;
                    else
                        asv1(i, 'genus') = cellstr('unclassified');
                        asv1(i, 'family') = cellstr('unclassified');
                        asv1(i, 'order') = cellstr('unclassified');
                        asv1(i, 'class') = cellstr('unclassified');
                        asv1(i, 'phylum') = cellstr('unclassified');
                    end
                end
            end
        end
    end
end

% now parse tables for each taxonomic level
taxLevels = asv1.Properties.VariableNames(5:10);
for i = 1:size(taxLevels,2)
    taxLevel = taxLevels{i};
    uniqueMembers = unique(asv1(:,taxLevel));
    subasv1 = asv1(:,{'Sample_ID', 'count', taxLevel});
    asv2 = unstack(subasv1, 'count', taxLevel);
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

asv2 = unstack(asv1(:,1:3),'count', 'ASV_asv');
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


