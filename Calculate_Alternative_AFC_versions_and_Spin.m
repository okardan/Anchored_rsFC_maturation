%########################
%# Title: Assessing neurocognitive maturation in early adolescence based on baby and adult functional brain landscapes
%# https://doi.org/10.1016/j.dcn.2025.101543
%# Contact: Omid Kardan omidk@med.umich.edu
%# This script calculates the (AFC) and also alternative versions of AFC in the posthoc analyses 

%# Requires parcel-wise preprocessed rsfMRI timeseries for ABCD participants at both Y0 and Y2 (for preprocessing details 
%# see fMRIprep supplementary material of this paper: Sripada, C., Angstadt, M., Taxali, A., Clark, D. A., Greathouse, T., Rutherford, S., ... & Heitzeg, M. 
%# 2021. Brain-wide functional connectivity patterns support general
%# cognitive ability and mediate effects of socioeconomic status in youth. Translational psychiatry, 11(1), 571.)

%# See also calculate_AFC_Kardan_et_al.m which is the script for the preprint version prior to peer reviews

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
%%
clear all
load('~\GAnetwork.mat')
load('~\GBnetwork.mat')
addpath('~')
folder_pish = '~\Rest\Gordon_Sub_Cere\Timeseries\';
T0 = readtable([folder_pish,'baselineYear1Arm1\ABCD_stats_postFD.csv']);
T2 = readtable([folder_pish,'2YearFollowUpYArm1\ABCD_stats_postFD.csv']);

outnet_edges_a = zeros(333);
num_nets_a = 13;

        
        for ii=1:333
            for jj=1:333
                    if GAnetwork(ii) ~= GAnetwork(jj)
                        outnet_edges_a(jj,ii) = 1;
                    end
            end
        end
   
innets_a = tril(ones(333),-1);
innets_a(outnet_edges_a == 1) =0;
 outnet_edges_a = tril(outnet_edges_a,-1);
 
 
outnet_edges_b = zeros(333);
num_nets_b = 11;

        
        for ii=1:333
            for jj=1:333
                if network(ii) ~= network(jj)
                    outnet_edges_b(jj,ii) = 1;
                end
            end
        end

innets_b = tril(ones(333),-1);
innets_b(outnet_edges_b == 1) =0;
 outnet_edges_b = tril(outnet_edges_b,-1);
idx = tril(ones(333),-1);

sum_r_vals =[]; subids =[]; Qs =[]; nruns =[]; pFD =[]; p_coefs = []; c_coefs =[];
names_blr = dir('~\connectomes\adult_net_arrange/*_baseline*rest*subc2024.mat');
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
        fname = ['~\connectomes\adult_net_arrange\',...
            char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 0 
        
        cnt = 0;
        mats = NaN(333,333,4); % 
        FDs = NaN(4,1);
        for run=1:4
            fname = ['~\connectomes\adult_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
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
        mats_b = NaN(333,333,4); % 
        for run=1:4
            fname = ['~\connectomes\baby_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_babywith_subc2024.mat'];
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
    
% original AFC
    sub_a_mat = squeeze(nanmean(mats,3));
    sub_b_mat = squeeze(nanmean(mats_b,3));
    sub_a_mat_mc =( atanh(sub_a_mat) - nanmean(atanh(sub_a_mat(find(idx)))) )/ std(atanh(sub_a_mat(find(idx))));
    sub_b_mat_mc =( atanh(sub_b_mat) - nanmean(atanh(sub_b_mat(find(idx)))) )/ std(atanh(sub_b_mat(find(idx))));

    
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_a==1))/5498 ...
        sum(sub_b_mat_mc(innets_b==1))/6017 ...
        sum(sub_a_mat_mc(outnet_edges_a==1))/49780 ...
        sum(sub_b_mat_mc(outnet_edges_b==1))/49261] ] ;

% modularity and AFC[modularity]
A = sub_a_mat(1:333,1:333); gamma = 1;
N=length(A);                            
K=sum(A);                               
m=sum(K);                               
B1=A-gamma*(K.'*K)/m;                   
Ci = GAnetwork;
s1=Ci(:,ones(1,N));                     
Qa1=~(s1-s1.').*B1/m;
Qa=sum(Qa1(:));


A = sub_b_mat(1:333,1:333); gamma = 1;
N=length(A);                            
K=sum(A);                               
m=sum(K);                               
B2=A-gamma*(K.'*K)/m;                    
Ci2 = network;
s2=Ci2(:,ones(1,N));                      
Qb1=~(s2-s2.').*B2/m;
Qb=sum(Qb1(:));
if ~isnan(Qb)    
[Ci,Qmax]=modularity_und(A,1);
numComm = length(unique(Ci));
else Qmax = NaN; numComm = NaN;
end

% AFC[participation coef]
A = sub_a_mat(1:333,1:333); Ci = GAnetwork;
n = length(A);
W_ = A.*(A>0);
S = sum(W_,2);
Gc = (W_~=0)*diag(Ci);
Sc2 = zeros(n,1);
for uu = 1:max(Ci)
    Sc2 = Sc2 + (sum(W_.*(Gc==uu),2).^2);
end
pc = ones(n,1) - Sc2./(S.^2);
pc(isnan(pc)) = 0;
pc(~pc) = 0;

A = sub_b_mat(1:333,1:333); Ci = network;
n = length(A);
W_ = A.*(A>0);
S = sum(W_,2);
Gc = (W_~=0)*diag(Ci);
Sc2 = zeros(n,1);
for uu = 1:max(Ci)
    Sc2 = Sc2 + (sum(W_.*(Gc==uu),2).^2);
end
pc_b = ones(n,1) - Sc2./(S.^2);
pc_b(isnan(pc_b)) = 0;
pc_b(~pc_b) = 0;

% AFC[clustering coefficient]
% local adult
cum_tul = 0; c_coefsA =[];
for loc = 1:max(GAnetwork)
    tul = length(find(GAnetwork == loc));
    Al = sub_a_mat(cum_tul+1:cum_tul+tul, cum_tul+1:cum_tul+tul);
    cum_tul = cum_tul + tul;
    W_ = Al.*(Al>0);
    K=sum(W_~=0,2);            	
cyc3=diag((W_.^(1/3))^3);           
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
c_coefsA = [c_coefsA; cyc3./(K.*(K-1))]; 
end
% local adult
cum_tul = 0; c_coefsB =[];
for loc = 1:max(network)
    tul = length(find(network == loc));
    Al = sub_b_mat(cum_tul+1:cum_tul+tul, cum_tul+1:cum_tul+tul);
    cum_tul = cum_tul + tul;
    W_ = Al.*(Al>0);
    K=sum(W_~=0,2);            	
cyc3=diag((W_.^(1/3))^3);           
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
c_coefsB = [c_coefsB; cyc3./(K.*(K-1))]; 
end



  nruns = [nruns; counter];
  pFD = [pFD; nanmean(FDs)];
    Qs = [Qs; [Qa Qb Qmax numComm]];
    c_coefs = [c_coefs; [nanmean(c_coefsA) nanmean(c_coefsB)  nanmean(c_coefG)]];
    p_coefs = [p_coefs; [nanmean(pc) nanmean(pc_b)]];
    subids = [subids; subid];
    k
end
wbdiff_adult = sum_r_vals(:,1)-sum_r_vals(:,3);
wbdiff_baby = sum_r_vals(:,2)-sum_r_vals(:,4);

T = table(subids,sum_r_vals,wbdiff_adult,wbdiff_baby,Qs,c_coefs, p_coefs, pFD, nruns);
writetable(T,['mod_sumrs_baseline_rest_mc2025_allruns.csv']);

% Y2
sum_r_vals =[]; subids =[]; Qs =[]; nruns=[]; pFD =[]; p_coefs =[]; c_coefs =[];
names_2yr = dir('~\connectomes\adult_net_arrange/*_2Year*rest*2024.mat');
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
        fname = ['~\connectomes\adult_net_arrange\',...
            char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
        if exist(fname)
            counter = counter+1;
        end
    end
    if counter > 0 
        
        cnt = 0;
        mats = NaN(333,333,4); % 
        FDs = NaN(4,1);
        for run=1:4
            fname = ['~\connectomes\adult_net_arrange\',...
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
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
            fname = ['~\connectomes\baby_net_arrange\',...
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_babywith_subc2024.mat'];
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
    
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_a==1))/5498 ...
        sum(sub_b_mat_mc(innets_b==1))/6017 ...
        sum(sub_a_mat_mc(outnet_edges_a==1))/49780 ...
        sum(sub_b_mat_mc(outnet_edges_b==1))/49261] ] ;
    
% modularity and AFC[modularity]
A = sub_a_mat(1:333,1:333); gamma = 1;
N=length(A);                            
K=sum(A);                               
m=sum(K);                               
B1=A-gamma*(K.'*K)/m;                    
Ci = GAnetwork;
s1=Ci(:,ones(1,N));                     
Qa1=~(s1-s1.').*B1/m;
Qa=sum(Qa1(:));


A = sub_b_mat(1:333,1:333); gamma = 1;
N=length(A);                           
K=sum(A);                           
m=sum(K);                              
B2=A-gamma*(K.'*K)/m;                   
Ci2 = network;
s2=Ci2(:,ones(1,N));                     
Qb1=~(s2-s2.').*B2/m;
Qb=sum(Qb1(:));
if ~isnan(Qb)    
[Ci,Qmax]=modularity_und(A,1);
numComm = length(unique(Ci));
else Qmax = NaN; numComm = NaN;
end

% AFC[participation coef]
A = sub_a_mat(1:333,1:333); Ci = GAnetwork;
n = length(A);
W_ = A.*(A>0);
S = sum(W_,2);
Gc = (W_~=0)*diag(Ci);
Sc2 = zeros(n,1);
for uu = 1:max(Ci)
    Sc2 = Sc2 + (sum(W_.*(Gc==uu),2).^2);
end
pc = ones(n,1) - Sc2./(S.^2);
pc(isnan(pc)) = 0;
pc(~pc) = 0;

A = sub_b_mat(1:333,1:333); Ci = network;
n = length(A);
W_ = A.*(A>0);
S = sum(W_,2);
Gc = (W_~=0)*diag(Ci);
Sc2 = zeros(n,1);
for uu = 1:max(Ci)
    Sc2 = Sc2 + (sum(W_.*(Gc==uu),2).^2);
end
pc_b = ones(n,1) - Sc2./(S.^2);
pc_b(isnan(pc_b)) = 0;
pc_b(~pc_b) = 0;

% AFC[clustering coefficient]
% local adult
cum_tul = 0; c_coefsA =[];
for loc = 1:max(GAnetwork)
    tul = length(find(GAnetwork == loc));
    Al = sub_a_mat(cum_tul+1:cum_tul+tul, cum_tul+1:cum_tul+tul);
    cum_tul = cum_tul + tul;
    W_ = Al.*(Al>0);
    K=sum(W_~=0,2);            	
cyc3=diag((W_.^(1/3))^3);           
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
c_coefsA = [c_coefsA; cyc3./(K.*(K-1))]; 
end
% local adult
cum_tul = 0; c_coefsB =[];
for loc = 1:max(network)
    tul = length(find(network == loc));
    Al = sub_b_mat(cum_tul+1:cum_tul+tul, cum_tul+1:cum_tul+tul);
    cum_tul = cum_tul + tul;
    W_ = Al.*(Al>0);
    K=sum(W_~=0,2);            	
cyc3=diag((W_.^(1/3))^3);           
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
c_coefsB = [c_coefsB; cyc3./(K.*(K-1))]; 
end



  nruns = [nruns; counter];
  pFD = [pFD; nanmean(FDs)];
    Qs = [Qs; [Qa Qb Qmax numComm]];
    c_coefs = [c_coefs; [nanmean(c_coefsA) nanmean(c_coefsB)  nanmean(c_coefG)]];
    p_coefs = [p_coefs; [nanmean(pc) nanmean(pc_b)]];
    subids = [subids; subid];
    k
end
wbdiff_adult = sum_r_vals(:,1)-sum_r_vals(:,3);
wbdiff_baby = sum_r_vals(:,2)-sum_r_vals(:,4);

T = table(subids,sum_r_vals,wbdiff_adult,wbdiff_baby,Qs,c_coefs, p_coefs, pFD, nruns);
writetable(T,['mod_sumrs_2Y_rest_mc2025_allruns.csv']);

%% Null models (Spin, Random, and Proximity)

% ###Spin### 
% requires the rotate_parcellation.m function from:
% Vasa F., Seidlitz J., Romero-Garcia R., Whitaker K. J., Rosenthal G., Ve?rtes P. E., Shinn M., 
% Alexander-Bloch A., Fonagy P., Dolan R. J., Goodyer I. M., the NSPN consortium, Sporns O., Bullmore E. T. (2017).
% Adolescent tuning of association cortex in human structural brain networks. Cerebral Cortex, 28(1):281–294.
% Updated in 2018 here: https://github.com/frantisekvasa/rotate_parcellation

clear all
addpath(genpath('Z:\abcd\UMich_connectomeDev\manuscript\DCN\Reviews\Spin_test'));
glr = importdata('Z:\abcd\UMich_connectomeDev\manuscript\DCN\Reviews\Spin_test/sphere_Gordon.txt');
coord_l = glr(1:161,:); coord_r = glr(162:333,:); nrot = 1000;
perm_id = rotate_parcellation(coord_l,coord_r,nrot);
% save perm_id
%
clear all
load('Z:\abcd\UMich_connectomeDev\manuscript\DCN\Reviews\Spin_test/perm_id.mat');
load('~\GAnetwork.mat')
load('~\GAmask.mat');
load('~\GBnetwork.mat')
load('~\GBmask.mat');
addpath('~')
folder_pish = '~\Rest\Gordon_Sub_Cere\Timeseries\';
Y02 = 'baseline';

for perm = 1:1000
    curr_perm = perm_id(:,perm);
ga_permnet = GAnetwork(curr_perm(GAmask));


    outnet_edges_a = zeros(333);
    num_nets_a = 13;
    
            for ii=1:333
                for jj=1:333
                    if ga_permnet(ii) ~= ga_permnet(jj)
                        outnet_edges_a(jj,ii) = 1;
                    end
                end
            end

    innets_a = tril(ones(333),-1);
    innets_a(outnet_edges_a == 1) =0;
    outnet_edges_a = tril(outnet_edges_a,-1);
    
    
gb_permnet = network(curr_perm(mask));
    outnet_edges_b = zeros(333);
    num_nets_b = 11;
    
            
            for ii=1:333
                for jj=1:333
                    if gb_permnet(ii)~= gb_permnet(jj)
                        outnet_edges_b(jj,ii) = 1;
                    end
                end
            end
 
    innets_b = tril(ones(333),-1);
    innets_b(outnet_edges_b == 1) =0;
    idx = tril(ones(333),-1);
    outnet_edges_b = tril(outnet_edges_b,-1);
    
    sum_r_vals =[]; subids =[]; 
    names_blr = dir(['~\connectomes\adult_net_arrange/*_',Y02,'*rest*subc2024.mat']);
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
            fname = ['~\connectomes\adult_net_arrange\',...
                char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
            if exist(fname)
                counter = counter+1;
            end
        end
        if counter > 0
            
                      cnt = 0;
            mats = NaN(333,333,4); %

            for run=1:4
                fname = ['~\connectomes\adult_net_arrange\',...
                    char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
                try
                    load(fname); adult_conn = full_conn;
                    cnt = cnt+1;
                    mats(:,:,run) = adult_conn(1:333,1:333);

                catch
                    mats(:,:,run) = NaN(333);

                    continue
                end
                
            end
            
            cnt = 0;
            mats_b = NaN(333,333,4); %
            for run=1:4
                fname = ['~\baby_net_arrange\',...
                    char(subid),'_baselineYear1Arm1_rest_',num2str(run),'_babywith_subc2024.mat'];
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
        
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_a==1))/5498 ...
            sum(sub_b_mat_mc(innets_b==1))/6017 ...
            sum(sub_a_mat_mc(outnet_edges_a==1))/49780 ...
            sum(sub_b_mat_mc(outnet_edges_b==1))/49261] ] ;
        

        subids = [subids; subid];
        k
    end
    wbdiff_adult = sum_r_vals(:,1)-sum_r_vals(:,3);
    wbdiff_baby = sum_r_vals(:,2)-sum_r_vals(:,4);
    
    T = table(subids,sum_r_vals,wbdiff_adult,wbdiff_baby);
    writetable(T,['~\Spin_test\rot_perms_gordonparcels/mod_sumrs_',Y02,'_rest_mc2024_',num2str(perm),'.csv']);
    
end

%
clear all
load('~\Spin_test/perm_id.mat');
load('~\GAnetwork.mat')
load('~\GAmask.mat');
load('~\GBnetwork.mat')
load('~\GBmask.mat');
addpath('~')
folder_pish = '~\Rest\Gordon_Sub_Cere\Timeseries\';
Y02 = '2YearFollowUpYArm1';

for perm = 1:100
    curr_perm = perm_id(:,perm);
ga_permnet = GAnetwork(curr_perm(GAmask));


    outnet_edges_a = zeros(333);
    num_nets_a = 13;
    
            for ii=1:333
                for jj=1:333
                    if ga_permnet(ii) ~= ga_permnet(jj)
                        outnet_edges_a(jj,ii) = 1;
                    end
                end
            end

    innets_a = tril(ones(333),-1);
    innets_a(outnet_edges_a == 1) =0;
    outnet_edges_a = tril(outnet_edges_a,-1);
    
    
gb_permnet = network(curr_perm(mask));
    outnet_edges_b = zeros(333);
    num_nets_b = 11;
    
            
            for ii=1:333
                for jj=1:333
                    if gb_permnet(ii)~= gb_permnet(jj)
                        outnet_edges_b(jj,ii) = 1;
                    end
                end
            end
 
    innets_b = tril(ones(333),-1);
    innets_b(outnet_edges_b == 1) =0;
    idx = tril(ones(333),-1);
    outnet_edges_b = tril(outnet_edges_b,-1);
    
    sum_r_vals =[]; subids =[]; 
    names_blr = dir(['~\adult_net_arrange/*_',Y02,'*rest*subc2024.mat']);
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
                char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
            if exist(fname)
                counter = counter+1;
            end
        end
        if counter > 0
            
            %    mats = NaN(333,333,4);
            cnt = 0;
            mats = NaN(333,333,4); %

            for run=1:4
                fname = ['~\adult_net_arrange\',...
                    char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_adultwith_subc2024.mat'];
                try
                    load(fname); adult_conn = full_conn;
                    cnt = cnt+1;
                    mats(:,:,run) = adult_conn(1:333,1:333);

                catch
                    mats(:,:,run) = NaN(333);

                    continue
                end
                
            end
            
            %    mats_b = NaN(333,333,4);
            cnt = 0;
            mats_b = NaN(333,333,4); %
            for run=1:4
                fname = ['~\baby_net_arrange\',...
                    char(subid),'_2YearFollowUpYArm1_rest_',num2str(run),'_babywith_subc2024.mat'];
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
        
        sum_r_vals = [sum_r_vals; [sum(sub_a_mat_mc(innets_a==1))/5498 ...
            sum(sub_b_mat_mc(innets_b==1))/6017 ...
            sum(sub_a_mat_mc(outnet_edges_a==1))/49780 ...
            sum(sub_b_mat_mc(outnet_edges_b==1))/49261] ] ;
        

        subids = [subids; subid];
        k
    end
    wbdiff_adult = sum_r_vals(:,1)-sum_r_vals(:,3);
    wbdiff_baby = sum_r_vals(:,2)-sum_r_vals(:,4);
    
    T = table(subids,sum_r_vals,wbdiff_adult,wbdiff_baby);
    writetable(T,['~\Spin_test\rot_perms_gordonparcels/mod_sumrs_',Y02,'_rest_mc2024_',num2str(perm),'.csv']);
    
end

% #####Null distributions from just adjacent voxels as clusters###
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


%% Plotting figures
clear all
T1 = importdata('~/mod_sumrs_randnull_baseline_rest_mc2024.csv');
T2 = importdata('~/mod_sumrs_randnull_2Y_rest_mc2024.csv');
T11 = readtable('~/mod_sumrs_anatnull_baseline_rest_mc2024.csv');
T22 = readtable('~/mod_sumrs_anatnull_2Y_rest_mc2024.csv');
abcdT0 = readtable('~\manuscript/mod_sumrs_baseline_rest_mc2024_allruns.csv');
abcdT2 = readtable('~\manuscript/mod_sumrs_2Y_rest_mc2024_allruns.csv');
abcdT1 = abcdT0(abcdT0.nruns > 1,:);

% Note: The Tu and Yeo networks are shared by the following studies:

% Tu, J. C., Myers, M., Li, W., Li, J., Wang, X., Dierker, D., ... & Wheelock, M. D. (2024). 
% Early life neuroimaging: The generalizability of cortical area parcellations across development. bioRxiv.

% Schaefer, Alexander, et al. Local-global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI.
% Cerebral cortex 28.9 (2018): 3095-3114.

%uiopen('~\Tu_Schaefer.csv',1);
% moo = TuSchaefer;
% moo0 = moo(string(moo.Session) == 'baselineYear1Arm1',:);
% moo2 = moo(string(moo.Session) == '2YearFollowUpYArm1',:);

WWs = [];
for i=1:1000
    ww = readtable(['~\Spin_test\rot_perms_gordonparcels/mod_sumrs_baseline_rest_mc2024_',num2str(i),'.csv']);
WWs = [WWs; [ww.wbdiff_adult ww.wbdiff_baby]];
i
end
    
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
h4 = histogram(WWs(:,1),5000)
h5 = histogram(WWs(:,2),5000)
h6 = histogram(T11.wbdiff_anatnull,50)
%h7 = histogram(moo0.diff_A,50)  % for plotting Yeo or Tu networks
%h8 = histogram(moo0.diff_B,50)  % for plotting Yeo or Tu networks
h1.FaceColor = [0.3 0.5 0.5]; h1.EdgeColor = [0.3 0.5 0.5];
h2.FaceColor = [0.3 0.4 0.75]; h2.EdgeColor = [0.3 0.4 0.75];
h3.FaceColor = [0.83 0.83 0.83]; h3.EdgeColor = [0.83 0.83 0.83];
h4.FaceColor = [0.5 0.5 0.5]; h4.EdgeColor = [0.5 0.5 0.5];
h5.FaceColor = [0.75 0.75 0.75]; h5.EdgeColor = [0.75 0.75 0.75];
h6.FaceColor = [0.83 0.68 0.68]; h6.EdgeColor = [0.83 0.68 0.68];
%h7.FaceColor = [0.9 0.6 0.73]; h7.EdgeColor = [0.9 0.6 0.73];
%h8.FaceColor = [0.6 0.9 0.68]; h8.EdgeColor = [0.6 0.9 0.68];
xlim([-0.2,1.84]);  title('abcd baseline');xlabel('Z_{within} - Z_{between}','FontSize',18); 
ylabel('Normalized Count','FontSize',18); ylim([0,1000]);
set(gca,'FontSize',18); legend('Adult Networks','Baby Networks','Null (Random)',...
    'Null (Adult Spin)','Null (Baby Spin)','Null (Proximity)',...%'Yeo17','Tu19',...
    'Location','northeast','NumColumns',2);

[a,b,c,stats] = ttest(abcdT1.wbdiff_baby,T11.wbdiff_anatnull) 
[a,b,c,stats] = ttest(abcdT1.wbdiff_baby,T11.wbdiff_anatnull) 
[a,b,c,stats] = ttest(abcdT2.wbdiff_adult,T22.wbdiff_anatnull) 
[a,b,c,stats] = ttest(abcdT2.wbdiff_adult,T22.wbdiff_anatnull) 


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


WWs2 = [];
for i=1:1000
    ww = readtable(['~\Spin_test\rot_perms_gordonparcels/mod_sumrs_2YearFollowUpYArm1_rest_mc2024_',num2str(i),'.csv']);
WWs2 = [WWs2; [ww.wbdiff_adult ww.wbdiff_baby]];
i
end
figure
subplot(3,1,1)
h1 = histogram(abcdT1.wbdiff_adult,50)
hold on
h2 = histogram(abcdT1.wbdiff_baby,50)
h3 = histogram(null_with_base - null_betw_base,5000)
h4 = histogram(WWs2(:,1),5000)
h5 = histogram(WWs2(:,2),5000)
h6 = histogram(T11.wbdiff_anatnull,50)
%h7 = histogram(moo2.diff_A,50)
%h8 = histogram(moo2.diff_B,50)
h1.FaceColor = [0.3 0.5 0.5]; h1.EdgeColor = [0.3 0.5 0.5];
h2.FaceColor = [0.3 0.4 0.75]; h2.EdgeColor = [0.3 0.4 0.75];
h3.FaceColor = [0.83 0.83 0.83]; h3.EdgeColor = [0.83 0.83 0.83];
h4.FaceColor = [0.5 0.5 0.5]; h4.EdgeColor = [0.5 0.5 0.5];
h5.FaceColor = [0.75 0.75 0.75]; h5.EdgeColor = [0.75 0.75 0.75];
h6.FaceColor = [0.83 0.68 0.68]; h6.EdgeColor = [0.83 0.68 0.68];
%h7.FaceColor = [0.9 0.6 0.73]; h7.EdgeColor = [0.9 0.6 0.73];
%h8.FaceColor = [0.6 0.9 0.68]; h8.EdgeColor = [0.6 0.9 0.68];
xlim([-0.2,1.84]);  title('abcd Y2');xlabel('Z_{within} - Z_{between}','FontSize',18); 
ylabel('Normalized Count','FontSize',18); ylim([0,1000]);
set(gca,'FontSize',18); legend('Adult Networks','Baby Networks','Null (Random)',...
    'Null (Adult Spin)','Null (Baby Spin)','Null (Proximity)',...%'Yeo17','Tu19',...
    'Location','northeast','NumColumns',2);
% export figure
