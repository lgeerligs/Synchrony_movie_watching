% This script reproduces part of the analyses in the paper “Age-related differences in 
% information processing during movie watching”. 
% https://doi.org/10.1016/j.neurobiolaging.2018.07.025

% In particular, it can be used to compute time-varying synchrony using sliding window analyses. 

% L Geerligs 14 sept 2018

clear all

datadir='/Data/';
resdir=[datadir '/results/'];

load([datadir 'Modules.mat'])
load([datadir 'Behavioral_data.mat'])
load([datadir 'residuals_subject_selection.mat'])
load([resdir 'synchrony_Craddock_smooth8_preprocWM_filterPW.mat'])

numTRs=192;
numSubjects=length(age);
winlengths=[15 20 30 45];
nWinlengths=length(winlengths);
nClusters=3;

%define covariates of no interest
covs_base=[rms_tot rms_max ICA_dat.rem];

%loop over window lengths
for ww=1:nWinlengths
    
    %get the sliding window
    windowlength=winlengths(ww);
    step=1;
    win=1:step:numTRs-windowlength;
    taper_weights=tukeywin(windowlength+2);
    taper_weights=taper_weights(2:end-1);
    
    %get timepoints in each window and their weight, so we can convert back
    %later
    timemat=zeros(numTRs, length(win));
    for w1 = 1:length(win)
        w=win(w1);
        timemat(w:w+windowlength-1,w1)=taper_weights;
    end
    
    %initialize variables
    dyn{1} = zeros(numSubjects, length(win),'single');
    dyn{2} = dyn{1} ;dyn{3} = dyn{1} ;
    
    %loop over the three clusters
    for net=1:nClusters
        disp([ww net])
        
        if net==1
            reg=FPNclust;
        elseif net==2
            reg=mPFCclust;
        elseif net==3
            reg=MTLclust;
        end
        
        %loop over subjects
        for i=1:numSubjects
          
            %this is the data of the participant
            subdata=squeeze(residuals(i,reg,:))';
            %get the average signal for all other participants
            j=setdiff(1:numSubjects,i);
            compdata=squeeze(mean(residuals(j,reg,:),1))';
            
            % iterate accross windows
            for w1 = 1:length(win)
                %get the current window
                w=win(w1);
                % get the current data
                tempdata=[subdata(w:w+windowlength-1,:) compdata(w:w+windowlength-1,:)];
                %compute the correlation
                r=weightedcorrs(tempdata,taper_weights');
                r=r(1:length(FPNclust),length(FPNclust)+1:end);
                r=rtoz(diag(r));
                %transform it with the Fisher r-z transformation
                dyn{net}(i,w1) =mean(r);
            end
        end

        %project synchorny estimates back to the original set of timepoints
        for i=1:numTRs
            vals=timemat(i,:);
            timedyn{net}(ww,:,i)=sum(vals'.*dyn{net}',1)./sum(vals);
        end
    end
end

%save results
save([resdir 'results_sliding_window_sync.mat'],'timedyn')



%% phase randomize and look timepoints with significant differences

nRandomizations=5000;
pval_thresh=0.025; % on each side (strong and weak age effects) 
names={'FPCN','mPFC','MTL'};

%initialize variables
mean_timedyn_pr_d=cell(nClusters,nWinlengths);
cage_timedyn_pr_s=mean_timedyn_pr_d; cage_timedyn=mean_timedyn_pr_d; page_timedyn=cage_timedyn;
times_highage=mean_timedyn_pr_d;times_lowage=mean_timedyn_pr_d;times_sig=mean_timedyn_pr_d;
netpairs_cor_mean=ones(nClusters,nWinlengths).*NaN;
pval_sim_netpairs_mean=netpairs_cor_mean;
pval_sim_netpairs_cage=netpairs_cor_mean;
netpairs_cor_cage=netpairs_cor_mean;

for ww=1:nWinlengths
    disp(ww)
    
    for net=1:nClusters
                
        %different phase randomization for different participants
        timedyn_pr_s=phaseran_diff(squeeze(timedyn{net}(ww,:,:))',nRandomizations);
        
        %same phase randomization for different participants
        timedyn_pr_d=phaseran(squeeze(timedyn{net}(ww,:,:))',nRandomizations);
        mean_timedyn_pr_d{net,ww}=squeeze(mean(timedyn_pr_d,2));
        
        for ran=1:nRandomizations
            %association between age and synchrony for the randomized data
            cage_timedyn_pr_s{net,ww}(ran,:)=partialcorr(squeeze(timedyn_pr_s(:,:,ran))', age,covs_base);
        end
        
        %association between age and synchrony for the real data
        [cage_timedyn{net,ww}, page_timedyn{net,ww}]=partialcorr(squeeze(timedyn{net}(ww,:,:)), age,covs_base);
      
        %timepoints with stronger or weaker than average age-effect
        times_highage{net,ww}=find(cage_timedyn{net,ww}>prctile(max(cage_timedyn_pr_s{net,ww}'),100-(pval_thresh.*100)));
        times_lowage{net,ww}=find(cage_timedyn{net,ww}<prctile(min(cage_timedyn_pr_s{net,ww}'),pval_thresh.*100));
        
        %timepoints with significant age-effect
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(page_timedyn{net,ww},0.05,'dep');
        times_sig{net,ww}=find(page_timedyn{net,ww}<=crit_p);

    end
    
    %look at the similarities betweeen the three networks
    for netpairs=1:nClusters
        if netpairs==1
            net1=1; net2=2;%FPCN & mPFC 
        elseif netpairs==2
            net1=1; net2=3;%FPCN & MTL
        elseif netpairs==3
            net1=2; net2=3;%mPFC & MTL
        end
        
        %similarity of the time-variation in mean-synchrony
        netpairs_cor_mean(netpairs,ww)=corr(squeeze(mean(squeeze(timedyn{net1}(ww,:,1:numTRs-1))',2)),squeeze(mean(squeeze(timedyn{net2}(ww,:,1:numTRs-1))',2)));
        ran_netpairs_cor_mean=diag(corr(mean_timedyn_pr_d{net1,ww}, mean_timedyn_pr_d{net2,ww}));
        pval_sim_netpairs_mean(netpairs,ww)=length(find(netpairs_cor_mean(netpairs,ww)<ran_netpairs_cor_mean))./nRandomizations;
        
        %similarity of the time-variation in age-effects
        netpairs_cor_cage(netpairs,ww)=corr(cage_timedyn{net1,ww}(1:numTRs-1), cage_timedyn{net2,ww}(1:numTRs-1));
        ran_netpairs_cor_cage=diag(corr(cage_timedyn_pr_s{net1,ww}', cage_timedyn_pr_s{net2,ww}'));
        pval_sim_netpairs_cage(netpairs,ww)=length(find(netpairs_cor_cage(netpairs,ww)<ran_netpairs_cor_cage))./nRandomizations;
    end
end


%save results
save([resdir 'results_sliding_window_sync.mat'],'cage_timedyn','page_timedyn','times_lowage','times_highage','times_sig','netpairs_cor_mean','pval_sim_netpairs_mean','netpairs_cor_cage','pval_sim_netpairs_cage', '-append')


%% plot results 
load([resdir 'results_sliding_window_sync.mat'])

% make plots of mean sync
for ww=2
    fig1=figure;hold on;
    shadedErrorBar(1:numTRs,mean(ztor(squeeze(timedyn{1}(ww,:,:)))),std(ztor(squeeze(timedyn{1}(ww,:,:))))./sqrt(numSubjects),'-r',1)
    shadedErrorBar(1:numTRs,mean(ztor(squeeze(timedyn{2}(ww,:,:)))),std(ztor(squeeze(timedyn{2}(ww,:,:))))./sqrt(numSubjects),'-b',1)
    shadedErrorBar(1:numTRs,mean(ztor(squeeze(timedyn{3}(ww,:,:)))),std(ztor(squeeze(timedyn{3}(ww,:,:))))./sqrt(numSubjects),'-c',1)
    xlabel('TR');ylabel('inter-subject synchrony (r)');
    fig1.Renderer='Painters';
    saveas(fig1,[resdir '/synchrony_all_' num2str(winlengths(ww)) '.pdf'])
end

%make plots of age-effects
colors={[1 0 0], [0 0 1], [0 1 0]};
for ww=2
    figure; hold on;
    for net=1:nClusters
        color=colors{net};
        for posneg=1:2
            if posneg==1
                alltimes=times_lowage{net,ww};
                alpha=0.15;
            elseif posneg==2
                alltimes=times_highage{net,ww};
                alpha=0.35;
            end
            if ~isempty(alltimes)
                sdec=diff(alltimes);
                point=sort([min(alltimes) alltimes(find(sdec>1)) alltimes(find(sdec>1)+1) max(alltimes)]);
                intervals=zeros(length(point)./2,2);
                vals=[1:2:length(point)];
                for i=1:length(point)./2
                    intervals(i,1)=point(vals(i));
                    intervals(i,2)=point(vals(i)+1);
                    
                    h=patch([intervals(i,1) intervals(i,2) intervals(i,2) intervals(i,1)], [-0.4 -0.4 0.2 0.2 ]',color);set(h,'EdgeColor','none','FaceAlpha',alpha)
                end
            end
        end
    end
    
    color={'r','b','c'};
    for net=1:nClusters
        solidline=cage_timedyn{net,ww}; dashedline=cage_timedyn{net,ww};
        nsig=setdiff(1:numTRs,times_sig{net,ww});
        solidline(nsig)=NaN;
        dashedline(times_sig{net,ww})=NaN;
        plot(dashedline,[':' color{net}],'LineWidth',2);
        plot(solidline,['-' color{net}],'LineWidth',2);
    end
    
    axis([1 numTRs -0.41 0.21]);ylabel('correlation synchrony age');xlabel('TR');set(gca,'Xtick',0:50:200);
    saveas(gcf,[resdir '/synchrony_ageeff_mPFCFPNMTL_' num2str(winlengths(ww)) '.pdf'])
end
