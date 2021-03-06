% This script reproduces part of the analyses in the paper �Age-related differences in 
% information processing during movie watching�. 
% https://doi.org/10.1016/j.neurobiolaging.2018.07.025

% In particular, it can be used to find communities of participants based on their cluster timecourses.

% L Geerligs 14 sept 2018

% this script requires the WSBM toolbox - http://tuvalu.santafe.edu/~aaronc/wsbm/  

clear all

datadir='/Data/';
resdir=[datadir '/results/'];

load([datadir 'Modules.mat'])
load([datadir 'residuals_subject_selection.mat'])
load([resdir 'synchrony_Craddock_smooth8_preprocWM_filterPW.mat'])

%combine data from all three clusters
meanTCs={meanTCs_mPFC', meanTCs_MTL',meanTCs_FPN'};
names={'mPFC_participant_clusters','MTL_participant_clusters','FPN_participant_clusters'};
sync={cluster_sync_mPFC,cluster_sync_MTL,cluster_sync_FPN};
names2={'mPFC','MTL','FPN'};

%% perform wsbm community detection

%loop over the three clusters
for net=1:3
    cmat=rtoz(corr(meanTCs{net}));
    cmat(eye(size(cmat))==1)=0;
    %loop over the number of communities that need to be detected
    for i=2:8
        [labels{i},model{i}] = wsbm(cmat,i,'alpha',0,'parallel',1,'numTrials',250);
    end
    
    save([resdir names{net}],'labels','model')
end

%% visualize results

age_communities_m={};age_communities_sd={};

plotvals=4;

%loop over the three clusters
for net=1:3
    load([resdir names{net} '.mat'])
    name=names2{net};
    cluster_sync=sync{net};
    cmat=rtoz(corr(meanTCs{net}));
    cmat(eye(size(cmat))==1)=0;
    sn=length(age);
    h=figure;
    
    %loop over the number of communities that need to be detected
    for i=2:8
        
        %get the mean synchrony in each community
        lab=labels{i};
        meansync=zeros(i,1);
        for j=1:i
            meansync(j)=mean(cluster_sync(lab==j));
        end
        
        %order communities by synchrony
        [~,sorder]=sort(meansync,'descend');
        label=lab;
        for j=1:i
            label(lab==sorder(j))=j;
        end
        [vals,sorder]=sort(label);
        
        %get the age per cluster
        for jj=1:i
            age_communities_m{i}(net,jj)=mean(age(label==jj));
            age_communities_sd{i}(net,jj)=std(age(label==jj));
        end
        
        %plot the correlation matrix of all values of k in one figure
        figure(h); cur=subplot(2,5,i-1);imagesc(ztor(cmat(sorder,sorder)));colormap(jet); caxis([-0.4 0.4]);
        set(cur,'dataAspectRatio',[1 1 1],'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
        subborder=find(diff(vals)==1);hold on;
        for ii=1:length(subborder);
            plot([subborder(ii) subborder(ii)],[1 sn],'-k');plot([1 sn],[subborder(ii) subborder(ii)],'-k');
        end
        saveas(gcf,[resdir name '_participant_cluster_FCmat_allsizes'  '.pdf'])
        
        if i==plotvals
            %plot a single correlation matrix for k=4
            figure; imagesc(ztor(cmat(sorder,sorder)));colorbar; colormap(jet); caxis([-0.3 0.3])
            set(gca,'dataAspectRatio',[1 1 1],'FontSize',14,'Xticklabel',[],'Yticklabel',[]);
            subborder=find(diff(vals)==1);hold on;
            for ii=1:length(subborder)
                plot([subborder(ii) subborder(ii)],[1 sn],'-k');plot([1 sn],[subborder(ii) subborder(ii)],'-k');
            end
            saveas(gcf,[resdir name '_participant_cluster_FCmat'  '.pdf'])
            
            
            %plot the signals of all participants from each cluster
            figure; imagesc(meanTCs{net}(:,sorder)');colormap(jet);caxis([-3.5 3.5])
            set(gca,'FontSize',14,'Xtick',50:50:150,'Yticklabel',[]);
            subborder=find(diff(vals)==1);hold on;
            for ii=1:length(subborder)
                plot([1 192],[subborder(ii) subborder(ii)],'-k');
            end
            colorbar; xlabel('TR'); ylabel('Participants');

            %compare participant ages per cluster
            [pval{net},anovatab{net},stats{net}] = anova1(age, label, 'off');
            
            saveas(gcf,[resdir name  '_allparticipant_cluster_timecourses'  '.pdf'])
        end
    end
end
