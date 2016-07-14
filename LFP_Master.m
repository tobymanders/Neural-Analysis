% Run titleRAW
clear all;
clc;
close all;

%% PARAMETERS
params.blstart = -1;
params.blend = -0.5;
params.wire = 'micro';
params.region = 'hippocampus'; %Region of interest
params.denoiseMethod = 'PCA'; %choose between 'PCA' and 'notch'
params.pre = 2; %seconds before event
params.post = 3; %seconds after event
params.movingWin = [.8 0.02];
params.Fs = 1e3; % desired frequency for LFP
params.pad = 2; %
params.fpass = [0 20]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 1; % average over trials
params.err = [2 0.01]; % population error bars
params.notch = 60;
params.cutoff = 40;
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
h = waitbar(0,'Sit tight. Plundering data...');

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
    LFP3 = zeros(size(LFP2));
    for c = 1:channelNum
        display(sprintf('Channel %d...',c))
        LFP3(c,:) = filtfilt(b,a,LFP2(c,:));
    end
else
    display('Skipping low-pass filter...');
    LFP3 = LFP2;
end
%% Plot data
% Plot raw, unfiltered LFP data
display('Generating plots...');
for j = 1:channelNum
    ax1(j) = subtightplot(channelNum,3,j*3-2);
    plot(LFP(j,floor(.5*end):floor(.5*end+10000)));
end
hold on;
title(ax1(1),'10s of Raw Data');

%Plot newly cleaned data
for k = 1:channelNum
    ax2(k) = subtightplot(channelNum,3,k*3-1);
    plot(LFP2(k,floor(.5*end):floor(.5*end+10000)));
end
title(ax2(1),'10s of Denoised Data');

%Plot newly lowpassed data
for l = 1:channelNum
    ax3(l) = subtightplot(channelNum,3,l*3);
    plot(LFP3(l,floor(.5*end):floor(.5*end+10000)));
end

linkaxes([ax1(1:end) ax2(1:end) ax3(1:end)],'xy');
title(ax3(1),'10s of low-pass-filtered Data');
hold off;

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
    cfig = figure;
    for ch = 1:channelNum
        
        %% construct LFP tensor.
        for tt = 1:nTrials
            LFPmat(:,tt,ch,aS) = LFP3(ch,round(trialStarts(tt)-params.pre)*params.Fs:round(trialStarts(tt)+params.post)*params.Fs);
        end
        
%         %%remove discharge trials
        rejectTrials = outliers(squeeze(range(LFPmat(:,:,ch,aS))));
%         LFPmat(:,rejectTrials,ch,aS) = NaN;
        
        
        
        %% time vectors
        tsec = linspace(-params.pre,params.post,size(LFPmat,1));
        
        % plotting LFP and spectrograms
        handleF = ch;
        brdr = 0.04;
        %figure(handleF)
        % first, ERPs
        %plotmultipleaxes(2,1,4,brdr,handleF)
        
        axlfp(ch) = subtightplot(channelNum,1,ch);
        hold on
        plot(tsec,nanmean(LFPmat(:,trialType==1,ch,aS),2),'color',col0)
        plot(tsec,nanmean(LFPmat(:,trialType==2,ch,aS),2),'color',col1a)
        plot(tsec,nanmean(LFPmat(:,trialType==3,ch,aS),2),'color',col1b)
        plot(tsec,nanmean(LFPmat(:,trialType==4,ch,aS),2),'color',col2)
        %hold off
        %xlim([-1 2])
        %set(gca,'linewidth',2,'fontsize',16)
        %title(deblank(names{ch}))
        
    end
    title(axlfp(1),sprintf('aligned on %s',alignName));
    linkaxes(axlfp(1:channelNum));
    display (sprintf('SAVING figure %d...', aS));
    saveas(cfig,strcat('C:\Users\melete2\Desktop\TobyData\7-8-16_LFP_figs\','_',upper(params.region),'_', nvName, alignName, upper(params.wire)), 'pdf');
    waitbar(.1*aS/2, h);
    
end
display ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Part 1 Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%% [20160610] calculating trial averaged spectrograms.
display('Calculating SPECTROGRAMS for...');

%% Calculate spectrograms

for aS = 1:2
    switch aS
        case 1
            alignName = 'Cue';
        case 2
            alignName = 'Response';
    end
    spectrograms(aS) = figure('Name',alignName);
    plotnum = channelNum*4;
    for  ch = 1:channelNum
        [blspec,t,f] = mtspecgramc(LFPmat(:,:,ch,1),params.movingWin,params);
        tspec = linspace(-params.pre,params.post,length(t));
        blsubspec = repmat(nanmean(blspec(tspec>params.blstart & tspec<params.blend,:)),length(tspec),1);
        display(sprintf('Channel %d...',ch));
        for conf = 1:4
            prog = ((ch+(conf-1)/4)*aS/channelNum/2)*.9 + .1;
            waitbar(prog,h)
            currplot = (ch-1)*4+conf;
            cols = 8;
            rows = ceil(plotnum/cols);
            [S{conf,ch,aS},t,f] = mtspecgramc(LFPmat(:,trialType==conf,ch,aS),params.movingWin,params);
            specax(currplot)=subtightplot(rows,cols,currplot);
            imagesc(tspec,f,(S{conf,ch,aS}./blsubspec)');
            xlim([-1 2])
            set(gca,'linewidth',2,'fontsize',12)
            axis xy square
            title(sprintf('%d-%d',ch,conf),'FontSize',5);
        end
    end
    set(specax(2:end), 'XTickLabel',[], 'YTickLabel',[]);
    set(specax(2:end), 'XTickLabel',[], 'YTickLabel',[]);

    axes = findall(spectrograms(aS),'type','axes');
    linkaxes (axes,'xy');
    hold off;
    
    display(sprintf('SAVING figure %d...',aS))
    saveas(spectrograms(aS),strcat(nevp(1:end-4),'LFP/',sprintf('CUBF_%s_%s%s_%s_%s_TM', nevp(end-6:end-5),...
        nvName(1:end-4), upper(params.region),upper(alignName),'_',upper(params.wire), date)),'pdf')

end %loop alignment
display ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Part 2 Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    close(h);
