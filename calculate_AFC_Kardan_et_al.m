%########################
%# Title: Assessing neurocognitive maturation in early adolescence based on baby and adult functional brain landscapes
%# Contact: Omid Kardan omidk@med.umich.edu
%# This script calculates the anchored rsFC maturation score (AFC) (also modularity for posthoc analyses)

%# Requires parcel-wise preprocessed rsfMRI timeseries for ABCD participants at both Y0 and Y2 (for preprocessing details 
%# see fMRIprep supplementary material of this paper: Sripada, C., Angstadt, M., Taxali, A., Clark, D. A., Greathouse, T., Rutherford, S., ... & Heitzeg, M. 
%# 2021. Brain-wide functional connectivity patterns support general
%# cognitive ability and mediate effects of socioeconomic status in youth. Translational psychiatry, 11(1), 571.)

%# The Supplemetary figures section requires the Matlab functions in the Parcel_Vis directory

clear all
folder_pish = '~\GordonParcelwise_timeseries\'; % location of rsfMRI parcelwise timeseries (1 file per run; up to 4 files per participant per year)

%%%%% 
% ABCD_rsfMRI_stats_FD.csv contains up to 4 rows per partitipant per year
% and inlcudes the motion and outlier information for each run from
% preprocessing step, as well as eventname (baseline or Y2)


T = readtable([folder_pish,'\ABCD_rsfMRI_stats_FD.csv']); num_rows = 25491+ 19199; % 25491 fMRI runs for Y0 and 19199 runs for Y2
load('~\GAmask.mat'); % for adult networks
load('~\GBmask.mat'); GBmask = mask;  % for baby networks

for k =1:num_rows
    if T.TRs(k) - T.confounds(k) < 250  % fewer than 250 DOF is excluded
        continue
    end
    subid = T.subjectkey{k};
    year = T.eventname{k};
    task = 'rest';
    run = T.run(k);
    try
        ptseries = importdata([folder_pish,year,'/Preprocessed/sub-',subid,'_ses-',year,'_task-',task,'_run-0',num2str(run),'_ts.txt']);
    catch
        continue
    end
    
    ptseries_gord = ptseries(:,1:333);
    adult_gord = ptseries_gord(:,GAmask);
    baby_gord = ptseries_gord(:,GBmask);
        
    full_conn = corr([adult_gord]);
    fname_a = ['~\adult_net_arrange\',subid,'_',year,'_',task,'_',num2str(run),'_adult.mat'];
    save(fname_a,'full_conn')
    
    
    full_conn = corr([baby_gord]);
    fname_b = ['~\baby_net_arrange\',subid,'_',year,'_',task,'_',num2str(run),'_baby.mat'];
    save(fname_b,'full_conn') 
 
end
%% calculating AFC 
clear all
load('~\GAnetwork.mat')
load('~\GBnetwork.mat')
addpath('E:\Omid\abcd\UMich_connectomeDev')
folder_pish = '~\GordonParcelwise_timeseries\';
T0 = readtable([folder_pish,'baselineYear1Arm1\ABCD_stats_postFD.csv']);
T2 = readtable([folder_pish,'2YearFollowUpYArm1\ABCD_stats_postFD.csv']);

% between and within net edges in adult arrangement
outnet_edges_a = zeros(333);
num_nets_a = 13;
for ntw1=1:num_nets_a-1
    for ntw2=ntw1+1:num_nets_a
        
        for ii=1:333
            for jj=1:333
                if GAnetwork(ii)==ntw1 & GAnetwork(jj)==ntw2
                    outnet_edges_a(jj,ii) = 1;
                end
            end
        end
    end
end
innets_a = tril(ones(333),-1);
innets_a(outnet_edges_a == 1) =0;

% between and within net edges in baby arrangement
outnet_edges_b = zeros(333);
num_nets_b = 11;
for ntw1=1:num_nets_b-1
    for ntw2=ntw1+1:num_nets_b
        
        for ii=1:333
            for jj=1:333
                if network(ii)==ntw1 & network(jj)==ntw2
                    outnet_edges_b(jj,ii) = 1;
                end
            end
        end
    end
end
innets_b = tril(ones(333),-1);
innets_b(outnet_edges_b == 1) =0;
idx = tril(ones(333),-1);

% read rearranged rsFC matrices from baseline year and calcualte [within -
% between]adult and [within - between]baby

sum_r_vals =[]; subids =[]; Qs =[]; nruns =[]; pFD =[];
names_blr = dir('~\adult_net_arrange/*_baseline*.mat');
subkeys =[];
for k =1:length(names_blr)
    subid1 = names_blr(k).name;
    subkey = subid1(1:15);
    subkeys = [subkeys; string(subkey)];
end
uq_subids = unique(subkeys);
for k =1:length(uq_subids)
subid = uq_subids(k);
    counter = 0;
    for run=1:4
        fname = ['~\adult_net_arrange\',...
            char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adult.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 0 
        cnt = 0;
        mats = NaN(333,333,4); 
        FDs = NaN(4,1);
        for run=1:4
            fname = ['~\adult_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adult.mat'];
            try
                load(fname); adult_conn = full_conn;
                 cnt = cnt+1;
                mats(:,:,run) = adult_conn(1:333,1:333);
                headm = T0.meanFD(string(T0.subjectkey) == subid & T0.run == run);
               FDs(run,1) = headm;
            catch
                mats(:,:,run) = NaN(333);
                FDs(run,1) = NaN;
                continue
            end

        end
        
        cnt = 0;
        mats_b = NaN(333,333,4); 
        for run=1:4
            fname = ['~\baby_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_baby.mat'];
            try
                load(fname); baby_conn = full_conn;
                cnt = cnt+1;
                mats_b(:,:,run) = baby_conn(1:333,1:333);
                
            catch
                mats_b(:,:,run) = NaN(333);
                continue
            end
            
        end
        
    end
    if counter ==0
        continue
    end
    sub_a_mat = squeeze(nanmean(mats,3));
    sub_b_mat = squeeze(nanmean(mats_b,3));
    sub_a_mat_mc =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));
    sub_b_mat_mc =( atanh(sub_b_mat) - nanmean(atanh(sub_b_mat(find(idx)))) )/ std(atanh(sub_b_mat(find(idx))));
% calculate mean connectivity within and mean connectivity between for baby and adult networks   
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_a==1))/5498 ...
        sum(sub_b_mat_mc(innets_b==1))/6017 ...
        sum(sub_a_mat_mc(outnet_edges_a==1))/49780 ...
        sum(sub_b_mat_mc(outnet_edges_b==1))/49261] ] ;
wbdiff_adult = sum_r_vals(:,1)-sum_r_vals(:,3);
wbdiff_baby = sum_r_vals(:,2)-sum_r_vals(:,4);
% calculate modularity (for posthoc analyses)
A = adult_conn(1:333,1:333); gamma = 1;
N=length(A);                            %number of vertices
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B1=A-gamma*(K.'*K)/m;                    %modularity matrix
Ci = GAnetwork;
s1=Ci(:,ones(1,N));                      %compute modularity
Qa1=~(s1-s1.').*B1/m;
Qa=sum(Qa1(:));
A = baby_conn(1:333,1:333); gamma = 1;
N=length(A);                            %number of vertices
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B2=A-gamma*(K.'*K)/m;                    %modularity matrix
Ci2 = network;
s2=Ci2(:,ones(1,N));                      %compute modularity
Qb1=~(s2-s2.').*B2/m;
Qb=sum(Qb1(:));
if ~isnan(Qb)    
[Ci,Qmax]=modularity_und(A,1);
numComm = length(unique(Ci));
else Qmax = NaN; numComm = NaN;
end
  nruns = [nruns; counter];
  pFD = [pFD; nanmean(FDs)];
    Qs = [Qs; [Qa Qb Qmax numComm]];
    subids = [subids; subid];
    k
end

T = table(subids,sum_r_vals,wbdiff_adult,wbdiff_baby,Qs,pFD, nruns);
writetable(T,['mod_sumrs_baseline_rest_mc2024_allruns.csv']); % This is used in the R script

% repeat for Y2
sum_r_vals =[]; subids =[]; Qs =[]; nruns=[]; pFD =[];
names_2yr = dir('~\adult_net_arrange/*_2Year*.mat');
subkeys =[];
for k =1:length(names_2yr)
    subid1 = names_2yr(k).name;
    subkey = subid1(1:15);
    subkeys = [subkeys; string(subkey)];
end
uq_subids = unique(subkeys);
for k =1:length(uq_subids)
subid = uq_subids(k);
    counter = 0;
    for run=1:4
        fname = ['~\adult_net_arrange\',...
            char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adult.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 0 
        
        cnt = 0;
        mats = NaN(333,333,4); 
        FDs = NaN(4,1);
        for run=1:4
            fname = ['~\adult_net_arrange\',...
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adult.mat'];
            try
                                load(fname); adult_conn = full_conn;
                 cnt = cnt+1;
                mats(:,:,run) = adult_conn(1:333,1:333);
                headm = T2.meanFD(string(T2.subjectkey) == subid & T2.run == run);
               FDs(run,1) = headm;
            catch
                mats(:,:,run) = NaN(333);
                FDs(run,1) = NaN;
                continue
            end
          
        end

        cnt = 0;
        mats_b = NaN(333,333,4); 
        for run=1:4
            fname = ['~\baby_net_arrange\',...
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_baby.mat'];
            try
                
                load(fname); baby_conn = full_conn;
                cnt = cnt+1;
                mats_b(:,:,run) = baby_conn(1:333,1:333);
                
            catch
                mats_b(:,:,run) = NaN(333);
                continue
            end
       
        end
        
    end
    if counter <=1
        continue
    end
      
    sub_a_mat = squeeze(nanmean(mats,3));
    sub_b_mat = squeeze(nanmean(mats_b,3));
    sub_a_mat_mc =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));
    sub_b_mat_mc =( atanh(sub_b_mat) - nanmean(atanh(sub_b_mat(find(idx)))) )/ std(atanh(sub_b_mat(find(idx))));
 % calculate mean connectivity within and mean connectivity between for baby and adult networks   
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_a==1))/5498 ...
        sum(sub_b_mat_mc(innets_b==1))/6017 ...
        sum(sub_a_mat_mc(outnet_edges_a==1))/49780 ...
        sum(sub_b_mat_mc(outnet_edges_b==1))/49261] ] ;
 wbdiff_adult = sum_r_vals(:,1)-sum_r_vals(:,3);
wbdiff_baby = sum_r_vals(:,2)-sum_r_vals(:,4);

% calculate modularity (for posthoc analyses)
    A = adult_conn(1:333,1:333); gamma = 1;
N=length(A);                            %number of vertices
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B1=A-gamma*(K.'*K)/m;                    %modularity matrix
Ci = GAnetwork;
s1=Ci(:,ones(1,N));                      %compute modularity
Qa1=~(s1-s1.').*B1/m;
Qa=sum(Qa1(:));
A = baby_conn(1:333,1:333); gamma = 1;
N=length(A);                            %number of vertices
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B2=A-gamma*(K.'*K)/m;                    %modularity matrix
Ci2 = network;
s2=Ci2(:,ones(1,N));                      %compute modularity
Qb1=~(s2-s2.').*B2/m;
Qb=sum(Qb1(:));
    
if ~isnan(Qb)    
[Ci,Qmax]=modularity_und(A,1);
numComm = length(unique(Ci));
else Qmax = NaN; numComm = NaN;
end
   pFD = [pFD; nanmean(FDs)];
   nruns = [nruns; counter];
    Qs = [Qs; [Qa Qb Qmax numComm]];
    
    subids = [subids; subid];
    k
end

T = table(subids,sum_r_vals,wbdiff_adult,wbdiff_baby, Qs,pFD, nruns);
writetable(T,['mod_sumrs_2Y_rest_mc2024_allruns.csv']);  % This is used in the R script


%% Null distributions from just adjacent voxels as clusters
clear all
cords1 = readtable('~/gordon_sub_cere_parcels_coords.csv');  % MNI coordinates of the Gordon parcels centers
cords = [cords1.Centroid_X(1:333), cords1.Centroid_Y(1:333), cords1.Centroid_Z(1:333)];
load('~\GAmask.mat');
[idx_cord,C,sumd,D] = kmeans(cords(GAmask,:),12);

innets_null = zeros(333); outnet_null = zeros(333);
num_nets_null = 12;  % This can also be set to 10 but results don;t change

for ntw1=1:num_nets_null
    innets_null(find(idx_cord == ntw1),find(idx_cord == ntw1)) = 1;
end
outnet_null(innets_null <1) = 1;
innets_null(find(tril(ones(333),-1)==0)) = 0;
outnet_null(find(tril(ones(333),-1)==0)) = 0;

idx = tril(ones(333),-1);

sum_r_vals =[]; subids =[]; Qs =[]; nruns =[];
names_blr = dir('~\adult_net_arrange/*_baseline*.mat');
subkeys =[];

for k =1:length(names_blr)
    subid1 = names_blr(k).name;
    subkey = subid1(1:15);
    subkeys = [subkeys; string(subkey)];
end
uq_subids = unique(subkeys);
for k =1:length(uq_subids)
subid = uq_subids(k);
    counter = 0;
    for run=1:4
        fname = ['~\adult_net_arrange\',...
            char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adult.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 1 
        
        mats = NaN(333,333,4);
        for run=1:4
            fname = ['~\adult_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adult.mat'];
            try
                load(fname); adult_conn = full_conn;
                mats(:,:,run) = adult_conn(1:333,1:333);
            catch
                mats(:,:,run) = NaN(333);
                continue
            end
        end
        
               
    end
    if counter <=1
        continue
    end
    sub_a_mat = squeeze(nanmean(mats,3)); 
    sub_a_mat_mc =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));


        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_null==1))/4725 ...
        sum(sub_a_mat_mc(outnet_null==1))/50553 ...
        ] ] ;
    
    subids = [subids; subid];
    k
end
wbdiff_anatnull = sum_r_vals(:,1)-sum_r_vals(:,2);

T = table(subids,sum_r_vals,wbdiff_anatnull);
writetable(T,['mod_sumrs_anatnull_baseline_rest_mc2024.csv']); % proximity-based (null) within and between networks mean connectivity

% repeat for Y2
sum_r_vals =[]; subids =[]; Qs =[]; nruns=[];
names_2yr = dir('~\adult_net_arrange/*_2Year*.mat');
subkeys =[];
for k =1:length(names_2yr)
    subid1 = names_2yr(k).name;
    subkey = subid1(1:15);
    subkeys = [subkeys; string(subkey)];
end
uq_subids = unique(subkeys);
for k =1:length(uq_subids)
subid = uq_subids(k);
    counter = 0;
    for run=1:4
        fname = ['\adult_net_arrange\',...
            char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adult.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 1 
        
        mats = NaN(333,333,4);
        for run=1:4
            fname = ['\adult_net_arrange\',...
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adult.mat'];
            try
                load(fname); adult_conn = full_conn;
                mats(:,:,run) = adult_conn(1:333,1:333);
            catch
                mats(:,:,run) = NaN(333);
                continue
            end
        end
        
               
    end
    if counter <=1
        continue
    end
    sub_a_mat = squeeze(nanmean(mats,3));
    
    sub_a_mat_mc =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));
        
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_null==1))/4725 ...
        sum(sub_a_mat_mc(outnet_null==1))/50553 ...
        ] ] ;
     
    subids = [subids; subid];
    k
end
wbdiff_anatnull = sum_r_vals(:,1)-sum_r_vals(:,2);

T = table(subids,sum_r_vals,wbdiff_anatnull);
writetable(T,['mod_sumrs_anatnull_2Y_rest_mc2024.csv']); % proximity-based (null) within and between networks mean connectivity

% random nuull dist
clear all
idx = tril(ones(333),-1);

sum_r_null_vals =[]; subids =[]; nruns =[];
names_blr = dir('~/*_baseline*.mat');
subkeys =[];
for k =1:length(names_blr)
    subid1 = names_blr(k).name;
    subkey = subid1(1:15);
    subkeys = [subkeys; string(subkey)];
end
uq_subids = unique(subkeys);
rng('shuffle');
for k =1:length(uq_subids)
subid = uq_subids(k);
    counter = 0;
    for run=1:4
        fname = ['~\adult_net_arrange\',...
            char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adult.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 1 
        
        mats = NaN(333,333,4);
        for run=1:4
            fname = ['~\adult_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adult.mat'];
            try
                load(fname); adult_conn = full_conn;
                mats(:,:,run) = adult_conn(1:333,1:333);
            catch
                mats(:,:,run) = NaN(333);
                continue
            end
        end
        
    end
    if counter <=1
        continue
    end
    sub_a_mat = squeeze(nanmean(mats,3));
    
    sub_a_mat_mc1 =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));
    sub_a_mat_mc = sub_a_mat_mc1(find(idx));
    
    sum_r_null_val1 =[];  sum_r_null_val2 =[];
    for j=1:100
        innets_null = randperm(55278,5498);
        outnet_null = setdiff(1:55278,innets_null);
        sum_r_null_val1 = [sum_r_null_val1  sum(sub_a_mat_mc(innets_null))/5498] ;
        sum_r_null_val2 = [sum_r_null_val2   sum(sub_a_mat_mc(outnet_null))/49780];
    end
 
    sum_r_null_vals = [sum_r_null_vals; sum_r_null_val1  sum_r_null_val2];
    subids = [subids; subid];
    k
end
with_null = sum_r_null_vals(:,1:100);
betw_null = sum_r_null_vals(:,201:200);

T = table(subids,with_null,betw_null);
writetable(T,['mod_sumrs_randnull_baseline_rest_mc2024.csv']);  % random networks at Y0


sum_r_null_vals =[]; subids =[]; 
names_2yr = dir('~\adult_net_arrange/*_2Year*.mat');
subkeys =[];
for k =1:length(names_2yr)
    subid1 = names_2yr(k).name;
    subkey = subid1(1:15);
    subkeys = [subkeys; string(subkey)];
end
uq_subids = unique(subkeys);
for k =1:length(uq_subids)
subid = uq_subids(k);
    counter = 0;
    for run=1:4
        fname = ['~\adult_net_arrange\',...
            char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adult.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 1 
        
        mats = NaN(333,333,4);
        for run=1:4
            fname = ['~\adult_net_arrange\',...
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adult.mat'];
            try
                load(fname); adult_conn = full_conn;
                mats(:,:,run) = adult_conn(1:333,1:333);
            catch
                mats(:,:,run) = NaN(333);
                continue
            end
        end
        
       
    end
    if counter <=1
        continue
    end
    sub_a_mat = squeeze(nanmean(mats,3));
       
    sub_a_mat_mc1 =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));
    sub_a_mat_mc = sub_a_mat_mc1(find(idx));
    
    sum_r_null_val1 =[];  sum_r_null_val2 =[];
    for j=1:100
        innets_null = randperm(55278,5498);
        outnet_null = setdiff(1:55278,innets_null);
        sum_r_null_val1 = [sum_r_null_val1  sum(sub_a_mat_mc(innets_null))/5498] ;
        sum_r_null_val2 = [sum_r_null_val2   sum(sub_a_mat_mc(outnet_null))/49780];
    end
 
    sum_r_null_vals = [sum_r_null_vals; sum_r_null_val1  sum_r_null_val2];
    subids = [subids; subid];
    k
end
with_null = sum_r_null_vals(:,1:100);
betw_null = sum_r_null_vals(:,101:200);

T = table(subids,with_null,betw_null);
writetable(T,['mod_sumrs_randnull_2Y_rest_mc2024.csv']); %random networks at Y2


%% Read the csv files generated above to make Figure 3
clear all
T1 = importdata('~/mod_sumrs_randnull_baseline_rest_mc2024.csv');
T2 = importdata('~/mod_sumrs_randnull_2Y_rest_mc2024.csv');
T11 = readtable('~/mod_sumrs_anatnull_baseline_rest_mc2024.csv');
T22 = readtable('~/mod_sumrs_anatnull_2Y_rest_mc2024.csv');
abcdT0 = readtable('~/mod_sumrs_baseline_rest_mc2024_allruns.csv');
abcdT2 = readtable('~/mod_sumrs_2Y_rest_mc2024_allruns.csv');
abcdT1 = abcdT0(abcdT0.nruns > 1,:);

null_with_base = reshape(T1.data(:,1:100),[6605*100,1]);
null_betw_base = reshape(T1.data(:,101:200),[6605*100,1]);
null_with_y2 = reshape(T2.data(:,1:100),[5170*100,1]);
null_betw_y2 = reshape(T2.data(:,101:200),[5170*100,1]);

figure
subplot(3,1,1)
h1 = histogram(abcdT1.wbdiff_adult,50)
hold on
h2 = histogram(abcdT1.wbdiff_baby,50)
h3 = histogram(null_with_base - null_betw_base,5000)
h4 = histogram(T11.wbdiff_anatnull,50)
h1.FaceColor = [0.3 0.5 0.5]; h1.EdgeColor = [0.3 0.5 0.5];
h2.FaceColor = [0.3 0.4 0.75]; h2.EdgeColor = [0.3 0.4 0.75];
h3.FaceColor = [0.73 0.73 0.73]; h3.EdgeColor = [0.73 0.73 0.73];
h4.FaceColor = [0.83 0.68 0.68]; h4.EdgeColor = [0.83 0.68 0.68];
xlim([-0.2,1.6]);  title('abcd baseline');xlabel('Conn_{within} - Conn_{between}','FontSize',18); 
ylabel('Normalized Count','FontSize',18); ylim([0,620]);
set(gca,'FontSize',18); legend('Adult Networks','Baby Networks','Null (Random)','Null (Proximity)');
addpath(genpath('~\export_fig-master')); % download export_fig from 
export_fig -m5 -transparent abcd_rest_bl_func_mat_null_hists2024.jpg

[a,b,c,stats] = ttest(abcdT1.wbdiff_baby,T11.wbdiff_anatnull) 
[a,b,c,stats] = ttest(abcdT1.wbdiff_baby,T11.wbdiff_anatnull) 
[a,b,c,stats] = ttest(abcdT2.wbdiff_adult,T22.wbdiff_anatnull) 
[a,b,c,stats] = ttest(abcdT2.wbdiff_adult,T22.wbdiff_anatnull) 

% Cohen's d values reported in Results section 3.1
numer = nanmean(abcdT1.wbdiff_baby) - nanmean(T11.wbdiff_anatnull);
denom = mean([nanstd(abcdT1.wbdiff_baby); nanstd(T11.wbdiff_anatnull)]);
cd_b = numer/denom  % 3.3867

numer = nanmean(abcdT1.wbdiff_adult) - nanmean(T11.wbdiff_anatnull);
denom = mean([nanstd(abcdT1.wbdiff_adult); nanstd(T11.wbdiff_anatnull)]);
cd_a = numer/denom  % 5.5705

numer = nanmean(abcdT2.wbdiff_baby) - nanmean(T22.wbdiff_anatnull);
denom = mean([nanstd(abcdT2.wbdiff_baby); nanstd(T22.wbdiff_anatnull)]);
cd_b = numer/denom  % 3.28

numer = nanmean(abcdT2.wbdiff_adult) - nanmean(T22.wbdiff_anatnull);
denom = mean([nanstd(abcdT2.wbdiff_adult); nanstd(T22.wbdiff_anatnull)]);
cd_a = numer/denom  % 5.79
 % Figure 3
figure
subplot(3,1,1)
h1 = histogram(abcdT2.wbdiff_adult,50)
hold on
h2 = histogram(abcdT2.wbdiff_baby,50)
h3 = histogram(null_with_y2 - null_betw_y2,5000)
h4 = histogram(T22.wbdiff_anatnull,50)
h1.FaceColor = [0.3 0.5 0.5]; h1.EdgeColor = [0.3 0.5 0.5];
h2.FaceColor = [0.3 0.4 0.75]; h2.EdgeColor = [0.3 0.4 0.75];
h3.FaceColor = [0.73 0.73 0.73]; h3.EdgeColor = [0.73 0.73 0.73];
h4.FaceColor = [0.83 0.68 0.68]; h4.EdgeColor = [0.83 0.68 0.68];
xlim([-.2,1.6]); title('abcd year2');xlabel('Conn_{within} - Conn_{between}','FontSize',18); 
ylabel('Normalized Count','FontSize',18); ylim([0,440]); 
set(gca,'FontSize',18);legend('Adult Networks','Baby Networks','Null (Random)','Null (Proximity)');
export_fig -m5 -transparent abcd_rest_y2_func_mat_null_hists2024.jpg
%% example sub connectome arrangements (Figure 2)
load('~\adult_net_arrange\NDARINVXXXXXX_baselineYear1Arm1_rest_2_adult.mat')
load('~\baby_net_arrange\NDARINVXXXXXX_baselineYear1Arm1_rest_2_baby.mat')
addpath(genpath('~\bluewhitered')) % download bluewhitered from 
figure;
subplot(1,2,1)
load('GBnetworknames.mat')
baby_conn1 = baby_conn;
baby_conn1(tril(ones(333),-1)==0) = 0;
imagesc(baby_conn1);colormap(bluewhitered), colorbar;
xnames =networknames(1:11); 
netsizes1 = [0 50 48 26 26 18 17 28 54 28 33];
netsizes = [50 48 26 26 18 17 28 54 28 33 5];
set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);

makans = cumsum(netsizes);
for j=1:length(netsizes)
line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
end
for j=1:length(netsizes)
line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
end
title('Average of rsFC (Baby network assignments)');
line([0 ,333],[0 ,333],'Color','black');
ylabel('average correlation value','Position',[381,167,1]); axis square


subplot(1,2,2)
adult_conn1 = adult_conn;
adult_conn1(tril(ones(333),-1)==0) = 0;
imagesc(adult_conn1);colormap(bluewhitered), colorbar;
load('GAnetnames.mat')
xnames = GAnetnames;
netsizes1 = [0 24    40     5    41    32    24    47     8    38     8     4    23];
netsizes = [24    40     5    41    32    24    47     8    38     8     4    23    39];
set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);

makans = cumsum(netsizes);
for j=1:length(netsizes)
line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
end
for j=1:length(netsizes)
line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
end

title('Average of rsFC (adult network assignments)');
line([0 ,333],[0 ,333],'Color','black');
ylabel('average correlation value','Position',[381,167,1]); axis square

%% Supplemetary figures (cortical surface maps)
% Script for Viewing Parcels on Surfaces
% The visualization functions are provided by Muriah Wheelock (mdwheelock@wustl.edu) written on 11-3-2020
clear all
addpath(genpath('~\Parcel_Vis'))

% Load your IM structure with desired network assignments %
load('~\Parcel_Vis\IM_11_BCP94_derivednetworks.mat') % load this for baby networks visualization
load('~\Parcel_Vis\IM_Gordon_2014_333_Parcels.mat') % load this for adult networks visualization
% Null distributions from just adjacent voxels as clusters
cords1 = readtable('~/gordon_sub_cere_parcels_coords.csv');
cords = [cords1.Centroid_X(1:333), cords1.Centroid_Y(1:333), cords1.Centroid_Z(1:333)];
load('~\GAmask.mat');
[idx_cord,C,sumd,D] = kmeans(cords(GAmask,:),12);
IM_null.name = 'anat_neighbor';
IM_null.key = [[1:333]' idx_cord];
IM_null.order = GAmask
IM_null.cMap = [0 1 1;.2 .6 .6;1 1 0;.5 0 0;0 0 .667;1 .7 .4;0 0.9 0;.69 .52 0;.8 .2 .8;1 0 0;.9 .9 .9;0 0 0];
IM_null.Nets = {'net1','net2','net3','net4','net5','net6','net7','net8','net9','net10','net11','net12'};
nullNets=IM_null.key(:,2);
nullorder = IM_null.order;
[B, I] = sort(nullorder);
nullNetOrigROIOrder = nullNets(I); % null networks ordered Left to Right

% now convert Parcels to Networks (in this example I used the null networks)
nullParcel_Nets.CtxL = zeros(32492,1);
nullParcel_Nets.CtxR = zeros(32492,1);
for i=1:333
    j = (Parcels.CtxL == i); %this is the Gordon vector indexing 333 order of nodes
    m = (Parcels.CtxR == i);
    for k=1:32492
        if j(k)==1
            nullParcel_Nets.CtxL(k)=nullNetOrigROIOrder(i);
        end
        if m(k)==1
        nullParcel_Nets.CtxR(k)=nullNetOrigROIOrder(i);
        end
    end
end

% View Parcels with coloring based on network
load('~\Parcel_Vis\Conte69_on_TT_32k.mat')
Anat.CtxL.data=nullParcel_Nets.CtxL; % plot desired networks
Anat.CtxR.data=nullParcel_Nets.CtxR;
% (2) Set parameters to view as desired
params.Cmap.P=IM_null.cMap;
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'
params.view='med';        % also, 'post','lat','med','dorsal'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);