% Run titleRAW


clear all;
clc;
close all;


%% PARAMETERS
params.testmode = 0;
params.band = 'all'; %choose b/w 'theta' and 'all'
params.blstart = -1;
params.blend = -0.5;
params.wire = 'micro';
params.region = 'hippocampus'; %Region of interest
params.denoiseMethod = 'notch'; %choose between 'PCA' and 'notch'
params.pre = 2; %seconds before event
params.post = 3; %seconds after event
params.movingWin = [.5 0.02];
params.Fs = 1e3; % desired frequency for LFP
params.pad = 2; %
params.fpass = [0 150]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 0; % average over trials
params.err = [2 0.01]; % population error bars
params.notch = 60;
params.cutoff = 200;
params.samplelen = 10; % in seconds
params.macronum = 2;
params.plotrange = [0 50];

% Define brain wave bands.
params.delta = [2 5];
params.theta = [5 9];
params.alpha = [9 13];
params.beta = [15 30];
params.lowgamma = [30 50];
params.highgamma = [75 150];
params.total = [2 150];

params.figdest = 'C:\Users\ShethLab\Desktop\Drive\FIGURES\';
switch params.wire
    case 'macro'
        params.source = '.ns3';
    case 'micro'
        params.source = '.ns5';
end

%% conflict colors
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;

%% Open NEV file
display ('Press ANY KEY to select your NEV file...');
pause;
[nevf, nevp] = uigetfile;
[~, nvName, ~] = fileparts([nevp nevf]);
NEV = openNEV([nevp nevf], 'nomat', 'nosave');
% h = waitbar(0,'Sit tight. Scrutinizing data...');

%% Pull out trial info
display ('Extracting trial information...');
% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);
% detect and remove false starts or practice
if max(diff(trigTimes) > 6)
    display(sprintf('\nIt seems there were practice trials in this file...\n\nRemoving %d events before event time gap.',find(diff(trigTimes) > 6)))
    trigs(1:find(diff(trigTimes) > 6)) = [];
    trigTimes(1:find(diff(trigTimes) > 6)) = [];
    nTrials = sum(trigs==90);
end

%% parsing behavior &  making a vector of conflict types.
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)

%% Parse other event times
respType = zeros(1,nTrials);
correct = trigs(trigs>=200 & trigs<=206);
respType(correct == 200 | correct == 204) = 1;    % CORRECT ANSWERS
respType(correct == 201 | correct == 205) = 2;   % INCORRECT ANSWERS


%% organizing responses & calculating reaction time.
responses = trigs(trigs>=100 & trigs<=104);
rt = trigTimes(trigs>=100 & trigs<104) - trigTimes(trigs>=1 & trigs<28);

%% Open NS5
display (sprintf('Opening raw %s data...', params.source));

if (exist([nevp nvName params.source]) == 2)
    NS5 = openNSx([nevp nvName params.source]);
else
    display(sprintf('No %s file found.', params.source));
end

%% Pull out, average and downsample LFP data for region of interest
% Find data for specified region
namesAll = {NS5.ElectrodesInfo.Label};
if strcmpi (params.region, 'hippocampus')
    regCh = find ((strncmp ('uHH',namesAll,3) == 1) | (strncmp ('uHT',namesAll,3) == 1)...
        | (strncmp ('uHB',namesAll,3) == 1) | (strncmp ('uHC',namesAll,3) == 1)| (strncmp ('uLHC',namesAll,3) == 1)...
        | (strncmp ('uRHC',namesAll,3) == 1) | (strncmp ('HH',namesAll,2) == 1) |(strncmp ('HT',namesAll,2) == 1)...
        |(strncmp ('HC',namesAll,2) == 1)|(strncmp ('LHC',namesAll,3) == 1));
elseif strcmpi (params.region, 'acc')
    regCh = find ((strncmp ('uAC',namesAll,3) == 1));
elseif strcmpi (params.region, 'all')
    regCh = find ((strncmp ('ainp',namesAll,4) ~= 1));
end
if strcmpi(params.wire, 'macro')
    regCh = regCh(params.macronum);
end
interval = NS5.MetaTags.SamplingFreq/params.Fs;
samples = length(NS5.Data);
names = namesAll(regCh);
channelNum = length(regCh);
% pull 1khz raw data from relevant channels
display ('Pulling data from NS5...');
LFP = zeros(channelNum,round(samples/interval));
if params.testmode
    channelNum = 1;
end
for i = 1:channelNum
    LFP(i,:) = NS5.Data(regCh(i),1:interval:end);
end
clear NS5;

%% Notch filter data
if strcmpi(params.wire, 'micro')
    display (sprintf('Removing noise using %s...\n',upper(params.denoiseMethod)));
    switch params.denoiseMethod
        case 'PCA'
            LFP2 = remove1stPC(LFP);
            display('...done.')
        case 'notch'
            display(sprintf('Applying NOTCH filter at %d Hz to...',params.notch));
            Wo = params.notch/(params.Fs/2);  BW = Wo/50;
            [b,a] = iirnotch(Wo,BW);
            LFP2 = zeros(size(LFP));
            for c = 1:channelNum
                display(sprintf('Channel %d...',c))
                LFP2(c,:) = filtfilt(b,a,LFP(c,:));
            end
        case 'chronotch'
            display('Applying CHRONUX NOTCH filter to...');
            LFP2 = zeros(size(LFP));
            for c = 1:channelNum
                display(sprintf('Channel %d...',c))
                LFP2(c,:) = rmlinesc(LFP(c,:)',params,[40 70])';
            end
        case 'matnotch'
            d = designfilt('bandstopiir','FilterOrder',4, ...
                'HalfPowerFrequency1',56,'HalfPowerFrequency2',64, ...
                'DesignMethod','butter','SampleRate',params.Fs);
            LFP2 = zeros(size(LFP));
            for c = 1:channelNum
                display(sprintf('Channel %d...',c))
                LFP2(c,:) = filtfilt(d,LFP(c,:));
            end
    end
else
    display('Skipped notch filter.')
    LFP2 = LFP;
end

%% Low pass filter data
if strcmp(params.wire,'micro')
    display('Applying LOW-PASS FILTER to...');
    Wn = params.cutoff/(0.5*params.Fs);
    [b,a] = butter (5, Wn);
    for c = 1:channelNum
        display(sprintf('Channel %d...',c))
        LFP3(c,:) = filtfilt(b,a,LFP2(c,:));
    end
else
    display('Skipped low-pass filter.');
    LFP3 = LFP2;
end

%% Band-pass for theta

LFP4 = zeros(size(LFP3));
%
% if strcmpi(params.band, 'theta')
%     display('Band-pass filtering for THETA on...');
%     Wn = params.theta/(params.Fs/2);
%     [z, p, k] = butter(40,Wn,'bandpass');
%     [sos,g] = zp2sos(z,p,k);
%     filt = dfilt.df2sos(sos,g);
%
%     for c = 1:channelNum
%         display(sprintf('Channel %d...',c))
%         LFP4(c,:) = filter(filt,LFP2(c,:));
%     end
% elseif strcmpi(params.band, 'all')
%         LFP4 = LFP3;
% else
%     error('Band not selected');
% end


LFP4 = zeros(size(LFP3));

if strcmpi(params.band, 'theta')
    display('Band-pass filtering for THETA on...');
    Wn = params.theta/(params.Fs/2);
    b = fir1(500,Wn,'bandpass');
    
    for c = 1:channelNum
        display(sprintf('Channel %d...',c))
        LFP4(c,:) = filtfilt(b,1,LFP2(c,:));
    end
elseif strcmpi(params.band, 'all')
    LFP4 = LFP3;
    display('Skipped band-pass for theta');
else
    error('Band not selected');
end


%% Generate sample plots

% Plot raw, unfiltered LFP data
display('Generating plots...');
fig(1) = figure;
for j = 1:channelNum
    ax1(j) = subtightplot(channelNum,3,j*3-2);
    plot(LFP(j,floor(.5*end):floor(.5*end+params.samplelen*params.Fs)));
end
hold on;
title(ax1(1),sprintf('%ds of Raw Data', params.samplelen));

%Plot newly cleaned data
for k = 1:channelNum
    ax2(k) = subtightplot(channelNum,3,k*3-1);
    plot(LFP2(k,floor(.5*end):floor(.5*end+params.samplelen*params.Fs)));
end
title(ax2(1),'10s of Denoised Data');

%Plot newly lowpassed data
for l = 1:channelNum
    ax3(l) = subtightplot(channelNum,3,l*3);
    plot(LFP3(l,floor(.5*end):floor(.5*end+params.samplelen*params.Fs)));
end

linkaxes([ax1(1:end) ax2(1:end) ax3(1:end)],'xy');
title(ax3(1),'10s of low-pass-filtered Data');
hold off;


%Plot theta
if strcmpi(params.band, 'theta')
    fig(2) = figure;
    for l = 1:channelNum
        ax4(l) = subtightplot(channelNum,1,l);
        plot(LFP4(l,floor(.5*end):floor(.5*end+10000)));
    end
    
    linkaxes(ax4(:),'xy');
    title(ax4(1),sprintf('%ds of theta',params.samplelen));
    hold off;
end

%% Plot data


% % Plot raw, unfiltered LFP data
%
% % alignName = 'Cue';
% % trialStarts =  trigTimes(trigs>=1 & trigs<28);
% % nTrials = sum(trigs>=1 & trigs<28);
% % display('Generating plots...');
% %
% % for j = 1:channelNum
% %     ax1(j) = subtightplot(channelNum,3,j*3-2);
% %     plot(LFP(j,1:floor(params.Fs*(trialStarts(1)+20))));
% %     for i=1:5
% %         line([trialStarts(i)*params.Fs trialStarts(i)*params.Fs], [-2000 2000], 'Color', 'Red');
% %     end
% % end
% % hold on;
% % title(ax1(1),'10s of Raw Data');
% %
% % %Plot newly cleaned data
% % for k = 1:channelNum
% %     ax2(k) = subtightplot(channelNum,3,k*3-1);
% %      plot(LFP2(j,1:floor(params.Fs*(trialStarts(1)+20))));
% %     for i=1:5
% %         line([trialStarts(i)*params.Fs trialStarts(i)*params.Fs], [-2000 2000], 'Color', 'Red');
% %     end
% % end
% % title(ax2(1),'10s of Denoised Data');
% %
% % %Plot newly lowpassed data
% % for l = 1:channelNum
% %     ax3(l) = subtightplot(channelNum,3,l*3);
% %      plot(LFP3(j,1:floor(params.Fs*(trialStarts(1)+20))));
% %     for i=1:5
% %         line([trialStarts(i)*params.Fs trialStarts(i)*params.Fs], [-2000 2000], 'Color', 'Red');
% %     end
% % end
% %
% % linkaxes([ax1(1:end) ax2(1:end) ax3(1:end)],'xy');
% % title(ax3(1),'10s of low-pass-filtered Data');
% % hold off;
% %




%% Plot spectrum for RT period vs. baseline
display('Generating TENSOR...');
LFPmat = zeros(((params.pre+params.post)*params.Fs)+1,nTrials,channelNum,2);

nTrials = sum(trigs>199 & trigs<207);
sub = ones(nTrials,1);

alignName{1} = 'Cue';
alignName{2} = 'Response';
alignName{3} = 'Feedback';
alignName{4} = 'Fixation';


for aS = 1:4
    switch aS
        case 1  %CUE
            trialStarts =  trigTimes(find((trigs>199 & trigs<207))-2*sub); 
        case 2  %RESPONSE
            trialStarts =  trigTimes(find((trigs>199 & trigs<207))-sub);
        case 3  %FEEDBACK
            trialStarts = trigTimes(trigs>199 & trigs<207);
        case 4  %FIXATION
            trialStarts = trigTimes(find((trigs>199 & trigs<207))-3*sub);
    end
    asfig(aS) = figure;
    rejectTrials = [];
    for ch = 1:channelNum
        %% construct LFP tensor.
        for tt = 1:nTrials
            LFPmat(:,tt,ch,aS) = LFP4(ch,round(trialStarts(tt)-params.pre)*params.Fs:round(trialStarts(tt)+params.post)*params.Fs);
        end
        
        %         %%remove discharge trials
        rejectTrials = [rejectTrials outliers(squeeze(range(LFPmat(:,:,ch,aS))))];
        
        
        
        
        %% time vectors
        tsec = linspace(-params.pre,params.post,size(LFPmat,1));
        
        % plotting LFP and spectrograms
        handleF = ch;
        brdr = 0.04;
        
        axlfp(ch) = subtightplot(channelNum,1,ch);
        hold on
        plot(tsec,nanmean(LFPmat(:,trialType==1,ch,aS),2),'color',col0)
        plot(tsec,nanmean(LFPmat(:,trialType==2,ch,aS),2),'color',col1a)
        plot(tsec,nanmean(LFPmat(:,trialType==3,ch,aS),2),'color',col1b)
        plot(tsec,nanmean(LFPmat(:,trialType==4,ch,aS),2),'color',col2)
        
    end
    
    % Delete the outlier channels
     uniqueRejects = unique(rejectTrials);
%     LFPmat(:,uniqueRejects,:,:) = NaN;
%     trialType(uniqueRejects)= NaN;
    
    title(axlfp(1),sprintf('aligned on %s BEFORE DISCHARGE REMOVAL',alignName{aS}));
    linkaxes(axlfp(1:channelNum));
    display (sprintf('SAVING figure %d...', aS));
    saveas(asfig(aS),strcat(params.figdest, sprintf('_%s_%s%s_%s_%s_', nevp(end-6:end-5),...
        nvName(1:end-4), upper(params.region),upper(alignName{aS}),upper(params.wire), date)),'pdf')
%     waitbar(.1*aS/2, h);
end
display (sprintf('DELETED %d out of %d trials due to discharges...',length(uniqueRejects),nTrials));
display ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Part 1 Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%% [20160610] calculating trial averaged spectrograms.
display('Calculating SPECTROGRAMS...');

%% Calculate spectrograms
% params.fpass = params.total;
bandname = 'total';

[~,t,f] = mtspecgramc(LFPmat(:,1,1,1),params.movingWin,params);
tspec = linspace(-params.pre,params.post,length(t));

parfor aS = 1:4
    for ch = 1:channelNum
        [nspec(aS,ch,:,:,:),~,~] = mtspecgramc(LFPmat(:,:,ch,aS),params.movingWin,params);
    end
end

%% Plot total spectrogram for each alignment

sp0 = squeeze(mean(nspec(1,:,:,:,:),2));
sp1 = mean(sp0,3);
bl1 = sp1(tspec>-1 & tspec<-.5,:);
bl2 = mean(bl1);
bl3 = repmat(bl2, size(sp1,1), 1);
bls = std(bl1);
bls1 = repmat(bls, size(sp1,1), 1);
figure;
for aS = 1:4
    sp0 = squeeze(mean(nspec(aS,:,:,:,:),2));
    sp1 = mean(sp0,3);
    zscore = (sp1 - bl3)./bls1;
    axh(aS) = subtightplot(1,4,aS);
    imagesc(tspec,f(f<50),zscore(:,f<50)', [-75 75]);
    axis xy square
    title(sprintf('%s',upper(alignName{aS})));
    colormap(jet)
end
linkaxes(axh(:))

%% Plot each trial 
figure;
nspect1 = squeeze(mean(nspec,1)) ;
rowcol = ceil(sqrt(length(nspect1(1,1,:))));

for trial = 1:100
    subtightplot(10,10,trial)
    imagesc(tspec,f,normlogspec(squeeze(nspect1(:,:,trial)))')
    title(trial)
    colormap(jet)
end

figure;
for trial = 101:200
    subtightplot(10,10,trial-100)
    imagesc(tspec,f,normlogspec(squeeze(nspect1(:,:,trial)))')
    title(trial)
    colormap(jet)
end

figure;
for trial = 201:nTrials
    subtightplot(10,10,trial-200)
    imagesc(tspec,f,normlogspec(squeeze(nspect1(:,:,trial)))')
    title(trial)
    colormap(jet)
end


%% Response Stuff
parfor ch = 1:channelNum
    [nsper(ch,:,:,:),~,~] = mtspecgramc(LFPmat(:,:,ch,2),params.movingWin,params);
end


spr1 = mean(nsper,1);
spr2 = mean(spr1,4);
spr3 = squeeze(spr2);

zscore = (spr3 - bl3)./bls1;
figure;
imagesc(tspec,f,zscore');
axis xy square

colormap(jet)

%% figs
for bS = 1:6
    specfig(bS) = figure;
end

%% Calculate spectrograms for ALL trials & channels and concatenate in third (trials) dimension



for aS = 1:2
    switch aS
        case 1
            alignName = 'Cue';
        case 2
            alignName = 'Response';
    end
    plotnum = channelNum*4;
    blspec = [];
    for  ch = 1:channelNum
        
        [blspec,t,f] = mtspecgramc(LFPmat(:,:,ch,1),params.movingWin,params);
        tspec = linspace(-params.pre,params.post,length(t));
        
        % Find loud baselines & delete those trials
        trange = tspec>params.blstart & tspec<params.blend;
        highbases = outliersright(squeeze(mean(mean(blspec(:,trange,:),2),1)));
        LFPmat(:,highbases,:,aS) = NaN;
        trialType(highbases) = NaN;
        blspec(:,:,highbases) = NaN;
        blspec2 = nanmean(blspec,3);
        
        blsubspec = repmat(mean(blspec2(tspec>params.blstart & tspec<params.blend,:)),length(tspec),1);
        cols = 2;
        rows = channelNum*4;
        display(sprintf('Channel %d... DONE!',ch));
        
        
        
        parfor conf = 1:4
            [S{conf,ch,aS},t,~] = mtspecgramc(LFPmat(:,trialType==conf,ch,aS),params.movingWin,params);
            U{conf,ch,aS} = squeeze(mean(S{conf,ch,aS},3));
            T(conf,ch,aS,:,:) = squeeze(mean(S{conf,ch,aS},3));
        end
        
        
        for bS = 1:6
            switch bS
                case 1
                    params.range = params.delta;
                    bandname = 'delta';
                case 2
                    params.range = params.theta;
                    bandname = 'theta';
                case 3
                    params.range = params.beta;
                    bandname = 'beta';
                case 4
                    params.range = params.lowgamma;
                    bandname = 'lowgamma';
                case 5
                    params.range = params.highgamma;
                    bandname = 'highgamma';
                case 6
                    params.range = params.alpha;
                    bandname = 'alpha';
            end
            set(specfig(bS),'name',bandname);
            set(0,'currentfigure',specfig(bS));
            hold on;
            fspec = linspace(params.fpass(1),params.fpass(2),length(blsubspec(1,:)));
            for conf = 1:4
                %             prog = ((ch+(conf-1)/4)*aS/channelNum/2)*.9 + .1;
                %             waitbar(prog,h)
                currplot = (ch-1)*4+conf;
                subplotnum = ((conf-1)*channelNum+ch)*2+aS-2;
                
                
                %             if ~strcmpi(params.band,'theta')
                %                 specax(currplot)=subtightplot(channelNum*2,4,currplot);
                %                 axis xy square
                %             else
                specax(currplot)=subtightplot(rows, cols, subplotnum);
                %             end
                imagesc(tspec,f,(U{conf,ch,aS}(:,fspec>params.range(1) & fspec<params.range(2))./blsubspec(:,fspec>params.range(1) & fspec<params.range(2)))');
                %imagesc(tspec,f,normlogspec((U{conf,ch,aS}))');
                
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',12)
                title(sprintf('%d-%d',ch,conf),'FontSize',5);
                
                %                 clear blspec
            end
            display (sprintf('%s band...',upper(bandname)));
            colormap(jet)

        end
        %     set(specax(2:end), 'XTickLabel',[], 'YTickLabel',[]);
        %         set(specax(1:end-1), 'XTickLabel',[], 'YTickLabel',[]);
        hold on;
    end %loop alignment
    
    
    
    %     axes = findall(fig(1),'type','axes');
    %     linkaxes (axes,'xy');
    hold off;
    display(sprintf('SAVING figure %d...',aS))
    
    
    
    %% Calculate power changes in different bands during task
    fpassrange = params.fpass(2)-params.fpass(1);
    [~,fbands,~] = size(S{1,1,1});
    
    
    
    low = floor((params.fpass(1)-params.fpass(1))/fpassrange*fbands+1);
    high = floor((params.fpass(2)-params.fpass(1))/fpassrange*fbands);
    
    fig(2) = figure;
    for ch = 1:channelNum
        for conf = 1:4
            trialpower = squeeze(nanmean(S{conf,ch,1}(:,low:high,:),2)); % power of theta in one trial over time
            trialnums = find(trialType==conf);
            for trial = 1:length(trialnums)
                rtpower(trial) = mean(trialpower(tspec>0 & tspec<rt(trialnums(trial)),trial)); % mean of theta power during rt
                basepower(trial) = mean(trialpower(tspec>params.blstart & tspec<params.blend,trial));
            end
            powerchange = 100*(rtpower-basepower)./basepower; % ratio of rt to bl power
            subplot(ceil(sqrt(channelNum*4)),ceil(sqrt(channelNum*4)),(ch-1)*4+conf)
            bar(sort(powerchange));
            title(sprintf('Channel %d, Conflict %d',ch,conf));
            meanchange(ch,conf) = nanmean(powerchange);
            semchange(ch,conf) = meanchange(ch,conf)/sqrt(length(powerchange));
            clear powerchange trialnums trialpower rtpower basepower
        end
    end
    
    
    fig(3) = figure;
    for ch = 1:channelNum
        subplot(2,4,ch);
        barwitherr(semchange(ch,:),meanchange(ch,:));
        title(sprintf('Channel %d',ch),'FontSize',8)
    end
    
    fig(4) = figure;
    barwitherr(mean(semchange,1),mean(meanchange,1));
    
    %Save figures
    parfor p = 1:4
        saveas(fig(p),strcat(params.figdest, sprintf('CUBF_%s_%s%s_%s_%s_', nevp(end-6:end-5),...
            nvName(1:end-4), upper(params.region),upper(alignName),upper(params.wire),upper(bandname), date),'_FIG_',num2str(p)),'pdf')
    end
end

%% Plot spectrogram for all trials, all conditions, all channels

figure;
allspec = squeeze(mean(squeeze(mean(T,1)),1));
bltot = mean(allspec(1,tspec>params.blstart & tspec>params.blend,:),2);
bltotrep = repmat(squeeze(bltot)',length(tspec),1);
imagesc(tspec,f,(squeeze(allspec(1,:,:))./bltotrep)'),
title('CUE ALL');
axis xy square
colormap(jet);

figure;
imagesc(tspec,f,(squeeze(allspec(2,:,:))./bltotrep)'),
title('RESPONSE ALL');
colormap(jet);
axis xy square


close(h);