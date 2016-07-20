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
params.movingWin = [.8 0.04];
params.Fs = 1e3; % desired frequency for LFP
params.pad = 2; %
params.fpass = [5 9]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 0; % average over trials
params.err = [2 0.01]; % population error bars
params.notch = 60;
params.cutoff = 50;
params.samplelen = 10; % in seconds

% Define brain wave bands.
params.delta = [2 5];
params.theta = [5 9];
params.alpha = [9 13];
params.beta = [15 30];
params.lowgamma = [30 50];
params.highgamma = [75 150];

params.figdest = 'C:\Users\melete2\Desktop\Drive\FIGURES\';
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
h = waitbar(0,'Sit tight. Scrutinizing data...');

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
        | (strncmp ('uRHC',namesAll,3) == 1) | (strncmp ('HH',namesAll,2) == 1) |(strncmp ('HT',namesAll,2) == 1));
elseif strcmpi (params.region, 'acc')
    regCh = find ((strncmp ('uAC',namesAll,3) == 1));
elseif strcmpi (params.region, 'all')
    regCh = find ((strncmp ('ainp',namesAll,4) ~= 1));
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


%% Low pass filter data
if strcmp(params.wire,'micro')
    display('Applying LOW-PASS FILTER to...');
    Wn = params.cutoff/(0.5*params.Fs);
    [b,a] = butter (2, Wn);
    for c = 1:channelNum
        display(sprintf('Channel %d...',c))
        LFP3(c,:) = filtfilt(b,a,LFP2(c,:));
    end
else
    display('Skipping low-pass filter...');
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


for aS = 1:2
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = sum(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
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
        uniqueRejects = unique(rejectTrials);
        LFPmat(:,uniqueRejects,:,aS) = NaN;
        trialType(uniqueRejects)= NaN;
        
        
        
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
    title(axlfp(1),sprintf('aligned on %s',alignName));
    linkaxes(axlfp(1:channelNum));
    display (sprintf('SAVING figure %d...', aS));
    saveas(asfig(aS),strcat(params.figdest, sprintf('_%s_%s%s_%s_%s_', nevp(end-6:end-5),...
            nvName(1:end-4), upper(params.region),upper(alignName),upper(params.wire), date)),'pdf')
    waitbar(.1*aS/2, h);
    
end
display ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Part 1 Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%% [20160610] calculating trial averaged spectrograms.
display('Calculating SPECTROGRAMS for...');

%% Calculate spectrograms

for bS = 1:5
    switch bS
        case 1
            params.fpass = params.delta;
            bandname = 'delta';
        case 2
            params.fpass = params.theta;
            bandname = 'theta';
        case 3
            params.fpass = params.beta;
            bandname = 'beta';
        case 4
            params.fpass = params.lowgamma;
            bandname = 'lowgamma';
        case 5
            params.fpass = params.highgamma;
            bandname = 'highgamma';
    end
    
    
    
    fig(1) = figure;
    
    
    
    for aS = 1:2
        switch aS
            case 1
                alignName = 'Cue';
            case 2
                alignName = 'Response';
        end
        plotnum = channelNum*4;
        for  ch = 1:channelNum
            [blspec,t,f] = mtspecgramc(LFPmat(:,~isnan(trialType),ch,1),params.movingWin,params);
            tspec = linspace(-params.pre,params.post,length(t));
            blspec2 = mean(blspec,3);
            blsubspec = repmat(mean(blspec2(tspec>params.blstart & tspec<params.blend,:)),length(tspec),1);
            display(sprintf('Channel %d...',ch));
            for conf = 1:4
%                 prog = ((ch+(conf-1)/4)*aS/channelNum/2)*.9 + .1;
                %             waitbar(prog,h)
                currplot = (ch-1)*4+conf;
                cols = 8;
                rows = ceil(plotnum/cols);
                [S{conf,ch,aS},t,f] = mtspecgramc(LFPmat(:,trialType==conf,ch,aS),params.movingWin,params);
                U{conf,ch,aS} = squeeze(mean(S{conf,ch,aS},3));
                %             if ~strcmpi(params.band,'theta')
                %                 specax(currplot)=subtightplot(channelNum*2,4,currplot);
                %                 axis xy square
                %             else
                specax(currplot)=subtightplot(channelNum*4,2,((conf-1)*channelNum+ch)*2+aS-2);
                %             end
                imagesc(tspec,f,(U{conf,ch,aS}./blsubspec)');
                %imagesc(tspec,f,normlogspec((U{conf,ch,aS}))');
                
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',12)
                title(sprintf('%d-%d',ch,conf),'FontSize',5);
                clear blspec
            end
        end
        %     set(specax(2:end), 'XTickLabel',[], 'YTickLabel',[]);
        set(specax(1:end-1), 'XTickLabel',[], 'YTickLabel',[]);
        hold on;
        
        
        
        %
    end %loop alignment
    colormap(jet)
    
    axes = findall(fig(1),'type','axes');
    linkaxes (axes,'xy');
    hold off;
    display(sprintf('SAVING figure %d...',aS))
    saveas(fig(1),strcat(params.figdest, sprintf('_%s_%s%s_%s_%s_', nevp(end-6:end-5),...
            nvName(1:end-4), upper(params.region),upper(alignName),upper(params.wire),upper(bandname), date)),'pdf')
    
    display (sprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %s Band Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',bandname));
 
    
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
    for p = 1:4
        saveas(fig(p),strcat(params.figdest, sprintf('CUBF_%s_%s%s_%s_%s_', nevp(end-6:end-5),...
            nvName(1:end-4), upper(params.region),upper(alignName),upper(params.wire),upper(bandname), date),'_FIG_',num2str(p)),'pdf')
    end
end
close(h);