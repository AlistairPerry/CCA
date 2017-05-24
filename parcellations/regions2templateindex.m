function [TempFinalIndex, Tempstrings] = regions2templateindex(parcniifile, tempniifile, tempLabels)

%regions2templateindex - For each parcellation region, find their best
%corresponding match within a target template image (i.e. Yeo 2011 Functional Clusters, AAL)

%[INPUT]:
    %parcniifile = filename of the source parcellation template
    %temppniifile = filename of the target template
        %!!! Both files must be in nii format and in the same space
    %tempLabels = strings corresponding to the integers within the target image 

%[OUTPUT]:
    %TempFinalIndex = Corresponding integer for each source region's best
    %match in the template image
    %TempFinalStrings = Corresponding label for the source region's best
    %match in the template image

%Alistair Perry(QIMR), 2017.

[parcdata] = load_untouch_nii(parcniifile);
numparcsints = max(max(max(parcdata.img)));

[tempdata] = load_untouch_nii(tempniifile);
numtempints = max(max(max(tempdata.img)));

for i = 1:numparcsints
    findparci = find(parcdata.img==i);
    tempmatch = tempdata.img(findparci);
    [tempmatches, ~, ~] = unique(tempmatch);
    tempmatches(tempmatches==0)=[];
    for j = 1:length(tempmatches);
        tempmatches(j,2) = length(find(tempmatch==tempmatches(j,1)));
    end
    if isempty(tempmatches)
        TempFinalIndex(i,1) = 0;
    else
        tempmatches = sortrows(tempmatches,2);
        TempFinalIndex(i,1) = tempmatches(j,1);
    end
end

Tempstrings = cell(numparcsints,1);
for i = 1:numtempints
    [Tempmatch,~] = find(TempFinalIndex==i);
    for j = 1:length(Tempmatch)
        Tempstrings(Tempmatch(j,1),1) = {[tempLabels{i,1} '_' int2str(j)]};
    end
end
end