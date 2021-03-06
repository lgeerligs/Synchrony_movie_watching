% This script reproduces part of the analyses in the paper �Age-related differences in 
% information processing during movie watching�. 
% https://doi.org/10.1016/j.neurobiolaging.2018.07.025

% In particular, it can be used to relate synchrony to functional connectivity. 
% It also computes the functional connectivity summary scores and relates those to the white matter (MK) values.  

% L Geerligs 14 sept 2018

% You will need to change the paths below to match your machine
% This script requires SPM12 - https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
% It also requires the Mediation Toolbox - https://github.com/canlab/MediationToolbox 
% as well as the Canlab core toolbox - https://github.com/canlab/CanlabCore 

clear all

datadir='/Data/';
resdir=[datadir '/results/'];

load([datadir 'Modules.mat'])
load([datadir 'Behavioral_data.mat'])
load([datadir 'residuals_subject_selection.mat'])
load([resdir 'synchrony_Craddock_smooth8_preprocWM_filterPW.mat'])

%define covariates of no interest
covs_base=[rms_tot rms_max ICA_dat.rem];

numSubjects=length(age);
nRegions=748;

%% compute functional connectivity

mat=ones(nRegions);
ind=find(triu(mat,1)==1);
FC=zeros(numSubjects, nRegions,nRegions);
for i=1:numSubjects
    disp(i)
    FC(i,:,:)=corr(squeeze(residuals(i,:,:))');
    meanFC(i)=mean(rtoz(FC(i,ind)));
end

%apply mean regression
FCmr=FC;
beta=[meanFC' ones(size(meanFC'))]\squeeze(FC(:,:));
FCmr(:,:)=squeeze(FC(:,:))-[meanFC' ones(size(meanFC'))]*beta;
FCmr(:,:)=squeeze(FCmr(:,:))+repmat(squeeze(nanmean(FC(:,:),1)),[length(meanFC) 1]);

%% compute stats and visualize results about synchrony relates to connectivity everywhere in the brain

%find the ROIs in each network which were not significantly related to age
BSnonage=setdiff(find(modcor==4), MTLclust);
FPNnonage=setdiff(find(modcor==14), FPNclust);
DMNnonage=setdiff(find(modcor==16), mPFCclust);

%define borders between modules for plotting
border=[];
border(1)=max(find(modcor==3));
border(2)=max(find(modcor==7));

%determine plot order and the borders for plotting
plotorder=[find(modcor<4)'; BSnonage'; MTLclust ;find(modcor>4&modcor<14)';  FPNclust; FPNnonage';  find(modcor==15)'; mPFCclust; DMNnonage'];
bordermPFC(1)=find(plotorder==min(mPFCclust));
bordermPFC(2)=find(plotorder==max(mPFCclust));
borderMTL(1)=find(plotorder==min(MTLclust));
borderMTL(2)=find(plotorder==max(MTLclust));
borderFPN(1)=find(plotorder==min(FPNclust));
borderFPN(2)=find(plotorder==max(FPNclust));

%mean connectivity
cmat1=squeeze(mean(FCmr));
cmat1(eye(size(cmat1))==1)=0;
figure; imagesc(cmat1(plotorder,plotorder));colorbar; colormap(jet); caxis([-1 1]);title('standard FC')
set(gca,'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
set(gca,'dataAspectRatio',[1 1 1]); hold on;
for i=1:length(border)
    plot([border(i) border(i)],[1 nRegions],'-k');plot([1 nRegions],[border(i) border(i)],'-k');
end
saveas(gcf,[resdir 'mean_standard_FC.pdf'])

%relate mPFC synchrony to connectivity
cmat=zeros(nRegions);pmat=cmat;
%correlate connectivity and synchrony
[cmat(:), pmat(:)]=partialcorr(FCmr(:,:), cluster_sync_mPFC,[cluster_sync_FPN cluster_sync_MTL age covs_base ],'rows','complete');
%find fdr-corrected significant correlations
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat(ind),0.05,'dep');
cmat(pmat>crit_p)=0;cmat(eye(size(cmat))==1)=0;
mPFC_cmat=cmat;

%relate FPN synchrony to connectivity
cmat=zeros(nRegions);pmat=cmat;
%correlate connectivity and synchrony
[cmat(:), pmat(:)]=partialcorr(FCmr(:,:), cluster_sync_FPN,[cluster_sync_mPFC cluster_sync_MTL age covs_base],'rows','complete');
%find fdr-corrected significant correlations
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat(ind),0.05,'dep');
cmat(pmat>crit_p)=0;cmat(eye(size(cmat))==1)=0;
FPN_cmat=cmat;

%relate MTL synchrony to connectivity
cmat=zeros(nRegions);pmat=cmat;
%correlate connectivity and synchrony
[cmat(:), pmat(:)]=partialcorr(FCmr(:,:), cluster_sync_MTL,[cluster_sync_FPN cluster_sync_mPFC age covs_base ],'rows','complete');
%find fdr-corrected significant correlations
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat(ind),0.05,'dep');
cmat(pmat>crit_p)=0;cmat(eye(size(cmat))==1)=0;
MTL_cmat=cmat;

%age effect on connectivity
cmat=zeros(nRegions);pmat=cmat;
%correlate connectivity and synchrony
[cmat(:), pmat(:)]=partialcorr(FCmr(:,:), age,[covs_base],'rows','complete');
%find fdr-corrected significant correlations
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat(ind),0.05,'dep');
cmat(pmat>crit_p)=0;cmat(eye(size(cmat))==1)=0;
age_cmat=cmat;

%plot results - relation mPFC connectivty and synchrony
fig1=figure; imagesc(mPFC_cmat(plotorder,plotorder));colorbar; colormap(jet); caxis([-0.4 0.4]);title('FC related to mPFC sync')
set(gca,'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
set(gca,'dataAspectRatio',[1 1 1]); hold on;
%plot borders
for i=1:length(border);
    plot([border(i) border(i)],[1 nRegions],'-k');plot([1 nRegions],[border(i) border(i)],'-k');
end
for i=1:length(bordermPFC);
    plot([bordermPFC(i) bordermPFC(i)],[1 nRegions],'--k');plot([1 nRegions],[bordermPFC(i) bordermPFC(i)],'--k');
end
myColorMap = jet(256);myColorMap(128:129,:) = 1;
colormap(myColorMap);colorbar
saveas(gcf,[resdir 'FC_related_to_mPFC_sync.pdf'])


%plot results - relation FPN connectivty and synchrony
figure; imagesc(FPN_cmat(plotorder,plotorder));colorbar; colormap(jet); caxis([-0.4 0.4]);title('FC related to FPN sync')
set(gca,'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
set(gca,'dataAspectRatio',[1 1 1]); hold on;
%plot borders
for i=1:length(border)
plot([border(i) border(i)],[1 nRegions],'-k');plot([1 nRegions],[border(i) border(i)],'-k');
end
for i=1:length(borderFPN);
    plot([borderFPN(i) borderFPN(i)],[1 nRegions],'--k');plot([1 nRegions],[borderFPN(i) borderFPN(i)],'--k');
end
myColorMap = jet(256);myColorMap(128:129,:) = 1;
colormap(myColorMap);
colorbar
saveas(gcf,[resdir 'FC_related_to_FPN_sync.pdf'])



%plot results - relation MTL connectivty and synchrony
figure; imagesc(MTL_cmat(plotorder,plotorder));colorbar; colormap(jet); caxis([-0.4 0.4]);title('FC related to MTL sync')
set(gca,'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
set(gca,'dataAspectRatio',[1 1 1]); hold on;
%plot borders
for i=1:length(border)
plot([border(i) border(i)],[1 nRegions],'-k');plot([1 nRegions],[border(i) border(i)],'-k');
end
for i=1:length(borderMTL);
    plot([borderMTL(i) borderMTL(i)],[1 nRegions],'--k');plot([1 nRegions],[borderMTL(i) borderMTL(i)],'--k');
end
myColorMap = jet(256);myColorMap(128:129,:) = 1;
colormap(myColorMap);
colorbar
saveas(gcf,[resdir 'FC_related_to_MTL_sync.pdf'])

%plot results - effect age on FC
fig1=figure; imagesc(age_cmat(plotorder,plotorder));colorbar; colormap(jet); caxis([-0.4 0.4]);title('age related to FC')
set(gca,'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
set(gca,'dataAspectRatio',[1 1 1]); hold on;
for i=1:length(border);
    plot([border(i) border(i)],[1 nRegions],'-k');plot([1 nRegions],[border(i) border(i)],'-k');
end
myColorMap = jet(256);myColorMap(128:129,:) = 1;
colormap(myColorMap);
colorbar
saveas(gcf,[resdir 'FC_related_to_age.pdf'])


%% construct a summary functional connectivty score using a leave-one-out approach

totconposFPN=zeros(numSubjects,1);totconnegFPN=totconposFPN;
totconposmPFC=zeros(numSubjects,1);totconnegmPFC=totconposmPFC;
totconposMTL=zeros(numSubjects,1);totconnegMTL=totconposmPFC;

if ~exist([resdir 'conscores.mat'])
    %use leave-one-out approach to get connectivity estimates for each
    %participant and for each of the three clusters
    parfor i=1:numSubjects
        disp(i)
        subset=setdiff(1:numSubjects,i);
        [cmat, pmat]=partialcorr(FCmr(subset,ind), cluster_sync_FPN(subset),[cluster_sync_mPFC(subset) cluster_sync_MTL(subset) covs_base(subset,:) age(subset)],'rows','complete');
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat,0.05,'dep');
        cmat(pmat>crit_p)=0;
        signconneg=find(cmat<0);
        signconpos=find(cmat>0);
        totconposFPN(i)=mean(FCmr(i,ind(signconpos)),2);
        totconnegFPN(i)=mean(FCmr(i,ind(signconneg)),2);
        
        [cmat, pmat]=partialcorr(FCmr(subset,ind), cluster_sync_mPFC(subset),[cluster_sync_FPN(subset) cluster_sync_MTL(subset) covs_base(subset,:) age(subset)],'rows','complete');
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat,0.05,'dep');
        cmat(pmat>crit_p)=0;
        signconneg=find(cmat<0);
        signconpos=find(cmat>0);
        totconposmPFC(i)=mean(FCmr(i,ind(signconpos)),2);
        totconnegmPFC(i)=mean(FCmr(i,ind(signconneg)),2);
        
        [cmat, pmat]=partialcorr(FCmr(subset,ind), cluster_sync_MTL(subset),[cluster_sync_FPN(subset) cluster_sync_mPFC(subset) covs_base(subset,:) age(subset)],'rows','complete');
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pmat,0.05,'dep');
        cmat(pmat>crit_p)=0;
        signconneg=find(cmat<0);
        signconpos=find(cmat>0);
        totconposMTL(i)=mean(FCmr(i,ind(signconpos)),2);
        totconnegMTL(i)=mean(FCmr(i,ind(signconneg)),2);
    end
    save([resdir 'conscores'],'totconposmPFC','totconnegmPFC','totconposFPN','totconnegFPN','totconposMTL','totconnegMTL')
end

%% run mediation analysis 

load([resdir 'conscores.mat'])
input={cluster_sync_FPN,cluster_sync_FPN,cluster_sync_mPFC,cluster_sync_mPFC,cluster_sync_MTL,cluster_sync_MTL};
mediator={totconposFPN,totconnegFPN,totconposmPFC,totconnegmPFC,totconposMTL,totconnegMTL};

for i=1:length(mediator)
    [paths, stats2] = mediation(zscore(age), zscore(input{i}), zscore(mediator{i}),'covs',[covs_base ]);
    table_paths(i,:)=paths;
    table_tvals(i,:)=stats2.t;
    table_pvals(i,:)=stats2.p;
    perc_mediation(i)=stats2.mean(5)./stats2.mean(4).*100;
end


%% relate functional connectivity to white matter data

load([datadir 'WMdata.mat'])
load([resdir 'conscores.mat'])
covs_WM=TIV';

%remove outliers and participants with stripes in diffusion data
WMdata_noout=WMdata;
WMdata_noout(stripeind>0.1,:)=NaN;
WMdata_noout(isnan(stripeind),:)=NaN;
for i=1:size(WMdata,2)
    in=find(~isnan(WMdata(:,i)));
    vals=WMdata(in,i);
    outlier=[find(vals>(mean(vals)+3*std(vals))) ;find(vals<(mean(vals)-3*std(vals)))];
    WMdata_noout(in(outlier),i)=NaN;
end

conscore={totconposmPFC,totconnegmPFC,totconposFPN,totconnegFPN,totconposMTL,totconnegMTL};
name={'posmPFC','negmPFC','posFPN','negFPN','posMTL','negMTL'};

%relate WM to connectivity
clear sync_WM_FC psync_WM_FC
for ii=1:length(conscore)
    for i=1:size(WMdata_noout,2)
        in=find(~isnan(WMdata_noout(:,i)));
        stat=regstats(conscore{ii}(in), zscore([WMdata_noout(in,i) zscore(WMdata_noout(in,i)).*zscore(age(in)) age(in) covs_base(in,:) covs_WM(in,:)  ]));
        sync_WM_FC(ii,i,:)=stat.tstat.t([2 3]);
        psync_WM_FC(ii,i,:)=stat.tstat.pval([2 3]);
        dfe(i)=stat.tstat.dfe;
    end
end
psync_WM_FC(psync_WM_FC>0.05./216)=1;


%plot and index effects
clear ry ro rm po pm py ry_spear py_spear po_spear ro_spear pm_spear rm_spear
for ii=1:length(conscore)
    [tracts,eff]=find(squeeze(psync_WM_FC(ii,:,:))<0.05);
    if ~isempty(tracts)
        for i=tracts
            in=find(~isnan(WMdata_noout(:,i)));
            young=find(age(in)<prctile(age(in),1/3*100));
            old=find(age(in)>prctile(age(in),2/3*100));
            middle=setdiff(1:length(in), [young; old]);
            
            ny(ii,i)=length(young);
            no(ii,i)=length(old);
            nm(ii,i)=length(middle);
            
            %get residuals for correlation
            ressyn=regstats(conscore{ii}(in), [covs_base(in,:) covs_WM(in,:) ]);
            ressyn=ressyn.r;
            resWM=regstats(WMdata_noout(in,i), [covs_base(in,:) covs_WM(in,:) ]);
            resWM=resWM.r;
            
            [ry(ii,i),py(ii,i)]=corr(ressyn(young), resWM(young));
            [ro(ii,i),po(ii,i)]=corr(ressyn(old), resWM(old));
            [rm(ii,i),pm(ii,i)]=corr(ressyn(middle), resWM(middle));
            
            [ry_spear(ii,i),py_spear(ii,i)]=corr(ressyn(young), resWM(young),'type','Spearman');
            [ro_spear(ii,i),po_spear(ii,i)]=corr(ressyn(old), resWM(old),'type','Spearman');
            [rm_spear(ii,i),pm_spear(ii,i)]=corr(ressyn(middle), resWM(middle),'type','Spearman');
        
            figure; scatter(ressyn(young), resWM(young),'filled'); lsline;
            hold on; scatter(ressyn(middle), resWM(middle), 'filled'); lsline;
            hold on; scatter(ressyn(old), resWM(old), 'filled'); l=lsline;
            p = findobj(gca);
            legend(p([5 8 10]),{'old','middle','young'},'Location','SouthEast')
            h=get(gca,'colororder');
            set(l(1),'LineWidth',2,'Color',h(3,:)); set(l(2),'LineWidth',2,'Color',h(2,:));set(l(3),'LineWidth',2,'Color',h(1,:))
            set(gca,'FontSize',25);
            ylabel(labels{i}); xlabel(name{ii}); axis([min(ressyn)+0.1*min(ressyn) max(ressyn)+0.1.*max(ressyn) min(resWM)+0.1.*min(resWM) max(resWM)+0.1.*max(resWM)])
            saveas(gcf,[resdir name{ii} '_' labels{i} '.pdf'])
        end
    end
end

  
%plot tracts on brain
hdr=spm_vol([datadir 'JHU-ICBM-tracts-maxprob-thr25-2mm.nii']);
tract_map=spm_read_vols(hdr);
for ii=1:length(conscore)
    tract_all=zeros(size(tract_map));
    [sigtracts,eff]=find(squeeze(psync_WM_FC(ii,:,:))<0.05);
    if ~isempty(sigtracts)
        for i=1:length(sigtracts)
            tract_all(tract_map==sigtracts(i))=sigtracts(i);
        end
        hdr.fname=[resdir name{ii} '_' labels{sigtracts(i)} '.nii'];
        spm_write_vol(hdr, tract_all);
    end
end
    
