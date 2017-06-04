function [CCAout] = cca_functional(connectivitymatrices, DM, motionFD, COG)

% Codes for CCA as performed in Perry et al., (2017, in review)
% Steps included: Normalization, functional decompostion, CCA and basic
% visualization output

%Input:
% connectivitymatrices - Functional matrices for all subjects
    %N.B: Matrix structure is in form i x i x s
        % where i is number of brain-regions, and s is number of subjects 
% DM - independent variates (i.e. non-imaging measures) used in CCA
    %.N.B: Matrix structure is in form s x k
        % where k is each non-imaging measure
% motionFD - mean framewise displacement (FD) of functional images
% COG - centroids for each parcellation region used in functional networks
    %.N.B: Matrix structure is in form i x 3

%Dependencies required in Matlab Path:
%FSLNets: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets
%PALM: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM

%Note: Codes have been implemented from Smith et al., (2015, Nat Neuroscience) 
    % http://www.fmrib.ox.ac.uk/datasets/HCP-CCA/
%& subsequentl Modified by Alistair Perry (QIMR Berghofer), 2017

%setup options
numtopcons=250; %number of connections most strongly expressed by CCA modes

%setup output structure
CCAout = struct; 

%save basic input to aid future analysis
CCAout.DM=DM;
CCAout.motionFD=motionFD;

%set up an identity matix for the purposes of extracting the functional connectomes upper triangle
parcnum=size(connectivitymatrices,1); %pull out resolution of parcellation template used

identmat=zeros(parcnum,parcnum);
niter=0;
for k = 1:parcnum
    for j = k+1:parcnum
        if k == parcnum; break; end;
        niter=niter+1;
        identmat(k,j)=niter;
    end
end

%reshape all subjects connectivity matrices, saving only upper triangle of
%connections
numsubjs=size(connectivitymatrices,3); % identify number of subjects

for i = 1:numsubjs
    niter=0;
    for k = 1:parcnum
        if k==parcnum end;
        for j = k+1:parcnum
            niter=niter+1;
            reshapeall(i,niter)=connectivitymatrices(k,j,i);
        end
    end
end

%Demeaning of matrices
CCAout.NET=reshapeall;
NET1=nets_demean(CCAout.NET);
NET1=NET1/std(NET1(:)); % no norm
amNET=abs(mean(CCAout.NET));
NET3=nets_demean(CCAout.NET./repmat(amNET,size(CCAout.NET,1),1));
NET3(:,amNET<0.1)=[]; %remove badly conditioned columns
NET3=NET3/std(NET3(:)); % norm by mean of columns
CCAout.grot=[NET3];

%Remove motion confounds
conf=palm_inormal(motionFD);    % Gaussianise
conf=nets_normalise([conf conf.^2]);
NETd=nets_demean(CCAout.grot-conf*(pinv(conf)*CCAout.grot));

%Extract principal eigenvectors of functional connectomes
%Note: Number of eigenvectors will depend on number of DM's in the analysis
eignum=length(CCAout.DM(:,1));
[CCAout.uu1,CCAout.ss1,vv1]=nets_svds(CCAout.NETd,eignum);

%Determine proportion variance explained by each eigenvector

for i = 1:eignum
CCAout.varexp(1,i)=[CCAout.ss1(i,i)*100./sum(sum(CCAout.ss1))];
end
CCAout.cumvarexp=cumsum(CCAout.varexp);

%Now the relative variance explained, when decomposing the functional edges
%into a less-reduced matrix

neweignum=floor(numsubjs*0.25);
[~,ss1,~]=nets_svds(CCAout.NETd,neweignum);
for i = 1:neweignum
CCAout.relvarexp(1,i)=[ss1(i,i)*100./sum(sum(ss1))];
end
CCAout.relcumvarexp=cumsum(CCAout.relvarexp);

%Perform Canonical Correlates Analysis
[CCAout.grotA, CCAout.grotB, CCAout.grotR, CCAout.grotU, CCAout.grotV, CCAout.grotstats]=canoncorr(CCAout.DM,CCAout.uu1);

%Loading of each subject measure onto each modes connectivity patterns
CCAout.conload=corr(CCAout.DM,CCAout.grotV);

%Expression of the original connectome edges within the CCA modes
%Note, some edges (those badly conditioned) have been removed
grotAAd = corr(CCAout.grotV(:,1),CCAout.grot)';
remove=find(amNET<0.1);
grotAAdfinal=zeros(1,length(amNET));
grotAAdfinal(1,remove(1,:))=1;
nonremove=find(amNET>=0.1);
grotAAd=grotAAd';

for i = 1:length(grotAAd)
    grotAAdfinal(1,nonremove(1,i))=grotAAd(1,i);
end
grotAAdfinal(grotAAdfinal==1)=0;

%Sort connections by their weighting
grotsort(1:length(grotAAdfinal),1)=1:length(grotAAdfinal);
grotsort(1:length(grotAAdfinal),2)=grotAAdfinal;
grotsort=sortrows(grotsort,2);

%Identify top positive connections
topgrotpos=grotsort(length(grotsort)-numtopcons+1:length(grotsort),:);
for i = 1:length(topgrotpos)
    [x,y]=find(identmat==topgrotpos(i,1));
    CCAout.toppostable{i,1}=x; %region i
    CCAout.toppostable{i,2}=y; %region j
    CCAout.toppostable{i,3}=topgrotpos(i,2); %correlation value
end

%Identify top negative connections
topgrotneg=grotsort(1:250,:);
for i = 1:length(topgrotpos)
    [x,y]=find(identmat==topgrotneg(i,1));
    CCAout.topnegtable{i,1}=x; %region i
    CCAout.topnegtable{i,2}=y; %region j
    CCAout.topnegtable{i,3}=topgrotneg(i,2); %correlation value
end

%Output toppositive and topnegative connections (i.e. top 250 edges in each direction) for Brain Net Viewer
%Will output only first mode for now - as mode significance can be performed
%through parametric and non-parametric options

%Write top positive connections
topgrotposmat=zeros(parcnum,parcnum);
for i = 1:length(topgrotpos)
    topgrotposmat(identmat==topgrotpos(i,1))=topgrotpos(i,2);
end
dlmwrite(['CCA' '_' int2str(numtopcons) 'topposcons' '_mode1' '.edge'], topgrotposmat, 'delimiter', '\t');

topgrotposmatsym = topgrotposmat;
for i = 1 : parcnum
    for j = i + 1 : parcnum
        if i == parcnum; break; end
        topgrotposmatsym(j,i)= topgrotposmat(i,j);
    end
end

%Calculate number of edges for each node in top positive connections (i.e. their degree)
topnodesdeg=degrees_und(topgrotposmatsym);
topnodes=unique(cell2mat(toppostable(:,1:2)));

for i = 1:parcnum
    search=find(topnodes==i);
    if ~isempty(search)
        nodesz(i,1)=topnodesdeg(1,i);
    else
        nodesz(i,1)=0;
    end
end

%Write nodal information of top negative connections for BNV
fid = fopen(['CCA' '_' 'nodes' int2str(numtopcons) 'topposcons' '_mode1' '.nodes'], 'wt');
for i = 1:parcnum
    fprintf(fid, '%f\t%f\t%f\t%d\t%d\t%s\n', COG(i,1), COG(i,2), COG(i,3), 1, nodesz(i,1), '~');
end

%Now top negative connections

%Write top negative connections
topgrotnegmat=zeros(parcnum,parcnum);
for i = 1:length(topgrotneg)
    topgrotnegmat(identmat==topgrotneg(i,1))=topgrotneg(i,2);
end
dlmwrite(['CCA' '_' int2str(numtopcons) 'topnegcons' '_mode1' '.edge'], topgrotnegmat, 'delimiter', '\t');

%Calculate number of edges for each node in top negative connections (i.e. their degree)
topgrotnegmatsym = topgrotnegmat;
for i = 1 : parcnum
    for j = i + 1 : parcnum
        if i == parcnum; break; end
        topgrotnegmatsym(j,i)= topgrotnegmat(i,j);
    end
end

topnodesdeg=degrees_und(topgrotnegmatsym);
topnodes=unique(cell2mat(topnegtable(:,1:2)));

for i = 1:parcnum
    search=find(topnodes==i);
    if ~isempty(search)
        nodesz(i,1)=topnodesdeg(1,i);
    else
        nodesz(i,1)=0;
    end
end

%Write nodal information of top negative connections for BNV
fid = fopen(['CCA' '_' 'nodes' int2str(numtopcons) 'topnegcons' '_mode1' '.nodes'], 'wt');
for i = 1:parcnum
    fprintf(fid, '%f\t%f\t%f\t%d\t%d\t%s\n', COG(i,1), COG(i,2), COG(i,3), 1, nodesz(i,1), '~');
end

end
