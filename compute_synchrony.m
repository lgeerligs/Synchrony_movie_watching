% This script reproduces part of the analyses in the paper �Age-related differences in 
% information processing during movie watching�. 
% https://doi.org/10.1016/j.neurobiolaging.2018.07.025

% In particular, it can be used to computes the inter-subject synchrony for each subject. 
% It also computes the correlation between synchrony and age, and performs the computation 
% of Bayes factors. It also includes a number of control analyses, as described in the 
% supplementary materials and the section �Sensory loss�. 

% L Geerligs 14 sept 2018

% You will need to change the paths below to match your machine
% This script requires SPM12 - https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 

clear all

datadir='/Data/';
resdir=[datadir '/results/'];

mkdir(resdir);
load([datadir 'residuals_subject_selection.mat'])
load([ datadir 'Modules.mat'])

%% compute inter-subject synchrony scores  (ISS)
% sync_groupavg is the default estimate of ISS
% sync_groupavg_beta and sync_matchavg are only used in the supplementary
% materials of the paper. 

nRegions=748;
nSubjects=length(age);
nTRs=192;

sync_groupavg=zeros(nSubjects, nRegions);corr_sub_groupavg_beta=sync_groupavg;
sync_matchavg=sync_groupavg;
parfor i=1:nSubjects
    disp(i)
    
    %this is the data of the participant
    subdata=squeeze(residuals(i,:,:));
    
    %get the average signal for all other participants
    j=setdiff(1:nSubjects,i);
    compdata=squeeze(mean(residuals(j,:,:),1));
    
    %compute the correlation and perform an r-z transformation
    sync_groupavg(i,:)=rtoz(diag(corr(subdata', compdata')));
    
    %compute synchrony based on a regession model
    for r=1:nRegions
        Bt=pinv([compdata(r,:)' ones(nTRs,1)])*subdata(r,:)';
        Btest(r)=Bt(1);
    end
    
    sync_groupavg_beta(i,:)=Btest;  Btest=[];
    
    %get the average signal for 20 age-matched participants
    [ages,order]=sort(abs(age-age(i)));
    group=order(1:21);
    j=setdiff(group,i);
    compdata=squeeze(mean(residuals(j,:,:),1));
    
    %compute the correlation and perform an r-z transformation
    sync_matchavg(i,:)=rtoz(diag(corr(subdata', compdata'))); 
end

%% compute statistics and make plots 

%define covariates of no interest
covs_base=[rms_tot rms_max ICA_dat.rem];

%relate synchrony and age
[cage_sync_groupavg, pcage_sync_groupavg]=partialcorr(sync_groupavg, age,covs_base);
regs_cage_sync_groupavg=find((pcage_sync_groupavg.*nRegions)<0.05);
thcage_sync_groupavg=cage_sync_groupavg;
thcage_sync_groupavg((pcage_sync_groupavg.*nRegions)>0.05)=0;
figure; plot(thcage_sync_groupavg);

%identify regions with significant mean synchrony
[h,p,ci,stats] = ttest(sync_groupavg);

%relate synchrony and age in age-matched analyses
[cage_sync_matchavg, pcage_sync_matchavg]=partialcorr(sync_matchavg, age,covs_base);
regs_cage_sync_matchavg=find((pcage_sync_matchavg.*nRegions)<0.05);
thcage_sync_matchavg=cage_sync_matchavg;
thcage_sync_matchavg((pcage_sync_matchavg.*nRegions)>0.05)=0;

%look at synchorny using beta values
[cage_sync_groupavg_beta, pcage_sync_groupavg_beta]=partialcorr(age, sync_groupavg_beta,covs_base);
regs_cage_sync_beta=find((pcage_sync_groupavg_beta.*nRegions)<0.05);
thcage_sync_groupavg_beta=cage_sync_groupavg_beta;
thcage_sync_groupavg_beta((pcage_sync_groupavg_beta.*nRegions)>0.05)=0;

%save results
save([resdir 'synchrony_Craddock_smooth8_preprocWM_filterPW'],'sync_groupavg','sync_matchavg', 'sync_groupavg_beta',.....
    'cage_sync_groupavg','pcage_sync_groupavg', 'thcage_sync_groupavg', 'regs_cage_sync_groupavg',.....
    'cage_sync_matchavg','pcage_sync_matchavg', 'thcage_sync_matchavg', 'regs_cage_sync_matchavg',.....
    'cage_sync_groupavg_beta','pcage_sync_groupavg_beta', 'thcage_sync_groupavg_beta', 'regs_cage_sync_beta')

%% link synchrony to RSFA (see supplementary materials)

RSFA=std(residuals,[],3);

%association between RSFA and synchrony
[r,p]=corr(RSFA, sync_groupavg);
r=diag(r);
p=diag(p);
sig=find(p<0.05./nRegions);
cRSFA_sync_groupavg=zeros(nRegions,1);
cRSFA_sync_groupavg(sig)=r(sig);

%association between RSFA and age
[r,p]=partialcorr(RSFA, age,covs_base);
sig=find(p<0.05./nRegions);
cRSFA_age=zeros(nRegions,1);
cRSFA_age(sig)=r(sig);

%association between synchrony and age, adjusting for RSFA
r=[]; p=[];
for i=1:nRegions
    [r(i), p(i)]=partialcorr(sync_groupavg(:,i), age,[covs_base RSFA(:,i)]);
end
thcage_sync_groupavg_aRSFA=r;
thcage_sync_groupavg_aRSFA((p.*nRegions)>0.05)=0;
figure; plot(thcage_sync_groupavg_aRSFA);


%% get bayes factors for association between age and synchrony in each brain region (see fig 1C)

%compute Bayes factors
for i=1:size(sync_groupavg,2)
    R2=regstats(zscore(sync_groupavg(:,i)),zscore([age covs_base]),'linear','rsquare');
    R1=regstats(zscore(sync_groupavg(:,i)),zscore(covs_base),'linear','rsquare');
    [bf10(i)] =jzs_partcorbf(sqrt(R1.rsquare),sqrt(R2.rsquare),3,4,nSubjects);
end


%define categories based on the Bayes factors
for i=1:nRegions
    if bf10(i)>100&&cage_sync_groupavg(i)<0
        Bayesevidence(i)=5;
    elseif bf10(i)>30&&cage_sync_groupavg(i)<0
        Bayesevidence(i)=4;
    elseif bf10(i)>10&&cage_sync_groupavg(i)<0
        Bayesevidence(i)=3;
    elseif bf10(i)>3&&cage_sync_groupavg(i)<0
        Bayesevidence(i)=2;
    elseif bf10(i)>1&&cage_sync_groupavg(i)<0
        Bayesevidence(i)=1;
    elseif bf10(i)<1&bf10(i)>(1/3)
        Bayesevidence(i)=-1;
    elseif bf10(i)<(1/3)&bf10(i)>(1/10)
        Bayesevidence(i)=-2;
    elseif bf10(i)<(1/10)&bf10(i)>(1/30)
        Bayesevidence(i)=-3;
    elseif bf10(i)<(1/30)&bf10(i)>(1/100)
        Bayesevidence(i)=-4;
    elseif bf10(i)<(1/100)
        Bayesevidence(i)=-5; 
    elseif bf10(i)>100&&cage_sync_groupavg(i)>0
        Bayesevidence(i)=10;
    elseif bf10(i)>30&&cage_sync_groupavg(i)>0
        Bayesevidence(i)=9;
    elseif bf10(i)>10&&cage_sync_groupavg(i)>0
        Bayesevidence(i)=8;
    elseif bf10(i)>3&&cage_sync_groupavg(i)>0
        Bayesevidence(i)=7;
    elseif bf10(i)>1&&cage_sync_groupavg(i)>0
        Bayesevidence(i)=6;
    end
end


%% define the three groups of regions for further analysis

%find the number of regions with significant age-differences in each module
mods_reg=[];
for i=1:16
    regs=find(modcor==i);
    mods_reg(i)=length(intersect(regs, regs_cage_sync_groupavg));
end

%get module labels for significant regions
sigmods=(thcage_sync_groupavg<0).*modcor';

%get the ROIs that are part the three networks that show the strongest age-effects
sigregsDMN=(thcage_sync_groupavg<0).*[1:nRegions]'.*double(modcor'==16);
sigregsFPN=(thcage_sync_groupavg<0).*[1:nRegions]'.*double(modcor'==14);
sigregsBrainstem=(thcage_sync_groupavg<0).*[1:nRegions]'.*double(modcor'==4);

%for the DMN, include only the ROIs within the mPFC
nonmPFC=[742 748 741 745 662 663 673];
mPFCclust=setdiff(sigregsDMN(sigregsDMN>0), nonmPFC);

%for the Brainstem, include only the ROIs within the MTL
nonMTL=[243 248];
MTLclust=setdiff(sigregsBrainstem(sigregsBrainstem>0), nonMTL);

%for the FPCN, include all ROIs
FPNclust=intersect(find(modcor==14), regs_cage_sync_groupavg);

%define visual and auditory clusters for the control analyses
visclust=intersect(find(modcor==1), regs_cage_sync_groupavg);
audclust=intersect(find(modcor==2), regs_cage_sync_groupavg);

%create map with the final clusters for visualizing on the brain
clustermap=zeros(size(sigmods));
clustermap(mPFCclust)=16;
clustermap(MTLclust)=4;
clustermap(FPNclust)=14;

% get mean timecourses per cluster
meanTCs_mPFC=squeeze(mean(residuals(:,mPFCclust,:),2));
meanTCs_FPN=squeeze(mean(residuals(:,FPNclust,:),2));
meanTCs_MTL=squeeze(mean(residuals(:,MTLclust,:),2));

% get the mean synchrony per cluster
cluster_sync_FPN=squeeze(mean(sync_groupavg(:,FPNclust),2));
cluster_sync_MTL=squeeze(mean(sync_groupavg(:,MTLclust),2));
cluster_sync_mPFC=squeeze(mean(sync_groupavg(:,mPFCclust),2));
cluster_sync_vis=squeeze(mean(sync_groupavg(:,visclust),2));
cluster_sync_aud=squeeze(mean(sync_groupavg(:,audclust),2));
cluster_sync_all=squeeze(mean(sync_groupavg(:,regs_cage_sync_groupavg),2));

save([resdir 'synchrony_Craddock_smooth8_preprocWM_filterPW'],'MTLclust','FPNclust','mPFCclust','visclust', 'audclust',.......
    'meanTCs_mPFC','meanTCs_FPN', 'meanTCs_MTL', 'cluster_sync_FPN', 'cluster_sync_MTL', 'cluster_sync_mPFC', 'cluster_sync_vis',.....
    'cluster_sync_aud', 'cluster_sync_all', '-append')

%% relate age and synchrony in each of the three clusters

%for each of the clusters get the correlation with age
[r, p]=partialcorr(cluster_sync_FPN, age,[covs_base])
[r, p]=partialcorr(cluster_sync_mPFC, age,[covs_base])
[r, p]=partialcorr(cluster_sync_MTL, age,[covs_base])

%for each of the clusters get the correlation with age, after adjusting for
%synchrony values in the auditory and visual clusters
[r, p]=partialcorr(cluster_sync_mPFC, age,[covs_base cluster_sync_vis cluster_sync_aud])
[r, p]=partialcorr(cluster_sync_FPN, age,[covs_base cluster_sync_vis cluster_sync_aud])
[r, p]=partialcorr(cluster_sync_MTL, age,[covs_base cluster_sync_vis cluster_sync_aud])

%for each of the clusters get the correlation with age, after adjusting for
%synchrony values in each of the signficant auditory and visual ROIs
[r, p]=partialcorr(cluster_sync_mPFC, age,[covs_base sync_groupavg(:,[visclust; audclust])])
[r, p]=partialcorr(cluster_sync_FPN, age,[covs_base sync_groupavg(:,[visclust; audclust])])
[r, p]=partialcorr(cluster_sync_MTL, age,[covs_base sync_groupavg(:,[visclust; audclust])])

%compute the correlation between the synchorony scores in the three clusters
[r, p]=corr(cluster_sync_FPN, cluster_sync_mPFC)
[r, p]=corr(cluster_sync_MTL, cluster_sync_FPN)
[r, p]=corr(cluster_sync_MTL, cluster_sync_mPFC)


%% create nifti images to visualize the results

%initalize the images
hdr=spm_vol([datadir '/Craddock_ROIs_included_ordered.nii']);
img=spm_read_vols(hdr);
img_ROIgroups=zeros(size(img));
ROIlist=unique(img);
ROIlist=ROIlist(ROIlist~=0);

for i=1:10
    if i==1; labels=thcage_sync_groupavg; name='age_effects_synchrony.nii';
    elseif i==2; labels=sigmods; name='age_effects_modules.nii';
    elseif i==3; labels=ztor(mean(sync_groupavg)); name='mean_sync.nii';
    elseif i==4; labels=thcage_sync_matchavg; name='agematch_effects_synchrony.nii';
    elseif i==5; labels=thcage_sync_groupavg_beta; name='agebeta_effects_synchrony.nii';
    elseif i==6; labels=Bayesevidence; name='Bayesevidence.nii';
    elseif i==7; labels=clustermap; name='three_clusters.nii';
    elseif i==8; labels=cRSFA_sync_groupavg; name='relation_RSFA_sync.nii';
    elseif i==9; labels=thcage_sync_groupavg_aRSFA; name='age_effects_synchrony_RSFAcorrected.nii';
    elseif i==10; labels=cRSFA_age; name='relation_RSFA_age.nii';
    end
    
    %for each ROI, find the corresponding voxels and insert the correct
    %statistic
    tel=0;
    for i=1:nRegions
        ind=find(img==ROIlist(i));
        tel=tel+1;
        img_ROIgroups(ind)=labels(tel);
    end
    
    %write the images
    hdr.dt=[16,0];
    hdr.private.dat.dtype='FLOAT32-LE';
    hdr.fname=[resdir name];
    spm_write_vol(hdr, img_ROIgroups);
    
end