% This script reproduces part of the analyses in the paper “Age-related differences in 
% information processing during movie watching”. 
% https://doi.org/10.1016/j.neurobiolaging.2018.07.025

% In particular, it can be used to relate synchrony to cognitive performance scores. 

% L Geerligs 14 sept 2018

% You will need to change the paths below to match your machine

clear all

datadir='/Data/';
resdir=[datadir '/results/'];

%load required data
load([datadir 'Behavioral_data.mat'])
load([datadir 'residuals_subject_selection.mat'])
load([resdir 'synchrony_Craddock_smooth8_preprocWM_filterPW.mat'])

%concatenate data from different clusters
synscore={cluster_sync_FPN,cluster_sync_mPFC,cluster_sync_MTL};
name={'FPCN','mPFC','MTL'};

%define covariates of no interest
covs_base=[rms_tot rms_max ICA_dat.rem];

%fieldnames of behavioral data 
fields=fieldnames(behdata);

%covary out synchrony scores from other clusters yes/no
covout=0; %if covout=1 synchrony scores from other cluster are covaried out below

%loop over clusters
for ii=1:length(synscore)
    ex=setdiff(1:length(synscore),ii);
    
    %loop over cognitive performance measures
    for i=1:length(fields)
        
        %find subjects for whom behavioral data is available
        in=find(~isnan(behdata.(fields{i})));
        
        %control for mRT/sdRT
        if strcmp(fields{i}, 'CSD') 
            excov=[zscore(behdata.('CRT')(in)) zscore(behdata.('CRT')(in)).*age(in)];
        else
            excov=[];
        end
        
        %covary out synchrony scores from other clusters yes/no
        if covout==1
            excov=[excov synscore{ex(1)}(in) synscore{ex(2)}(in)] ;
        end
        
        %build regression model
        stat=regstats(synscore{ii}(in), zscore([behdata.(fields{i})(in) zscore(behdata.(fields{i})(in)).*zscore(age(in))  age(in) covs_base(in,:) excov ]));
        
        %save results
        sync_beh(ii,i,:)=stat.tstat.t([2 3]);
        psync_beh(ii,i,:)=stat.tstat.pval([2 3]);
        model_r2(ii,i)=stat.rsquare;
        model_f(ii,i)=stat.fstat.f;
        model_pval(ii,i)=stat.fstat.pval;
        df(i)=stat.tstat.dfe;
    end
end

%check which are significant after Bonferonni correction
psync_beh(psync_beh>(0.05./numel(psync_beh)))=1;


% for the interaction effects, look within different subgroups and make a scatterplot of the results
clear r_subgroup p_subgroup
for ii=1:length(synscore)

    [a,b]=find(squeeze(psync_beh(ii,:,:))<0.5);
    if ~isempty(a)
        for i=a'

            in1=find(~isnan(behdata.(fields{i})));
            
            %define subgroups
            old=find(age(in1)>prctile(age(in1),2/3*100));
            young=find(age(in1)<prctile(age(in1),1/3*100));
            middle=setdiff(1:length(in1), [young; old]);
        
            %degrees of freedom
            df_subgroup{ii}(i,3)=length(old);
            df_subgroup{ii}(i,2)=length(middle);
            df_subgroup{ii}(i,1)=length(young);
            
            %get residuals 
            ressyn=regstats(synscore{ii}(in1), [covs_base(in1,:)]);
            ressyn=ressyn.r;
            resbeh=regstats(zscore(behdata.(fields{i})(in1)), [covs_base(in1,:) ]);
            resbeh=resbeh.r;

            %correlations
            [r_subgroup{ii}(i,3),p_subgroup{ii}(i,3)]=corr(ressyn(old), resbeh(old));
            [r_subgroup{ii}(i,1),p_subgroup{ii}(i,1)]=corr(ressyn(young), resbeh(young)) ;
            [r_subgroup{ii}(i,2),p_subgroup{ii}(i,2)]=corr(ressyn(middle), resbeh(middle)) ;
            
            %plot results
            figure; scatter(ressyn(young), resbeh(young),'filled'); lsline;
            hold on; scatter(ressyn(middle), resbeh(middle), 'filled'); lsline;
            hold on; scatter(ressyn(old), resbeh(old), 'filled'); l=lsline;

            p = findobj(gca);
            legend(p([5 8 10]),{'old','middle','young'},'Location','SouthEast')
            h=get(gca,'colororder');
            set(l(1),'LineWidth',2,'Color',h(3,:)); set(l(2),'LineWidth',2,'Color',h(2,:));set(l(3),'LineWidth',2,'Color',h(1,:))
            set(gca,'FontSize',25);
            ylabel(fields{i}); xlabel(name{ii}); axis([min(ressyn)-0.1*min(ressyn) max(ressyn)+0.1.*max(ressyn) min(resbeh)-0.1*min(resbeh) max(resbeh)+0.1.*max(resbeh)])
            saveas(gcf,[resdir name{ii} '_' fields{i} '.pdf'])
        end
    end
end


