function [coherenceStats] = analyzeMSITcoherence(patientID,sessionNum,nevFile)
%ANALYZEMSITCOHERENCE does spike-field coherence analysis on MSIT data.
%
%   [coherenceStats] = analyzeMSITcoherence(patientID,sessionNum,nevFile)
%   calculates spike-field coherograms for MSIT data in nevFile and its
%   associated ns(3) file. analyzeMSITcoherence saves statistics, and
%   plots results.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)


%% loading data from NEV file
display('loading action potential and local field potential data...')
[dataPath, nvName, nvExt] = fileparts(nevFile);

% defining nsFile.
nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);

% parsing files and loading data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
    if ~exist(nsFile,'file')
        [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
        NS3 = openNSx(fullfile(nsPath,nsFile));
    else
        NS3 = openNSx(nsFile);
    end
elseif strcmp(nvExt,'.mat')
    load(nevFile);
end


%% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% [20160608] ruling out practice sessions or false starts.
if max(diff(trigTimes) > 6)
    %     figure(1123)
    %     maximize(1123)
    %     hold on
    %     plot(diff(trigTimes))
    %     title('visualized trial structure via differences between event times')
    %     text(0,10,'<--practice trials')
    %     text(find(diff(trigTimes) > 6),20,'test trials-->')
    %     hold off
    
    display(sprintf('\nIt seems there were practice trials in this file...\n\nRemoving %d events before event time gap.',find(diff(trigTimes) > 6)))
    
    % adjusting number of events.
    %     find(diff(trigTimes) > 6);% hard-code justification: the max ITI is 5.017
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


%% timing (seconds)
pre = 2;
post = 3;


%% parsing channel labels
macroLabels = deblank({NS3.ElectrodesInfo.Label})';
numBFs = length(macroLabels)./8;
% BFlabels =

for bf = 1:numBFs
    tmp = char(deblank(macroLabels((bf)*8,:)));
    BFlabels{bf} = tmp(1:end-1);
end
BFlabels


%% defining Chronux parameters.
movingWin = [0.8 0.010];
params.Fs = 2e3; % sampling frequency for LFP
params.pad = 2; %
params.fpass = [0 200]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 0; % average over trials {CHANGES BELOW}
params.err = [2 0.01]; % population error bars


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


%% [20160610] denoising LFP
denoiseMethod = 'notch'
switch denoiseMethod
    case 'PCA'
        display('denoising using PCA...')
        dData = remove1stPC(double(NS3.Data));
        display('...done.')
    case 'notch'
        Wo = 60/(params.Fs/2);  BW = Wo/50;
        [b,a] = iirnotch(Wo,BW);
        %        freqz(b,a);
        for c = 1:size(NS3.Data,1)
            display(sprintf('applying notch filter to channel %d',c))
            dData(c,:) = filtfilt(b,a,double(NS3.Data(c,:)));
        end
end


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% initializing LFPmat
LFPlabels = {NS3.ElectrodesInfo.Label};
LFPmat = zeros(((pre+post)*params.Fs)+1,nTrials,length(LFPlabels),2);
% for loop to save multiple epochs
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
            trigTimes(trigs>=100 & trigs<104)
    end
    
    for ch = 1:length(LFPlabels)
        for tt = 1:nTrials
            
            %%~~~~~~~~~~~~~~~~~~~FROM OLD (WORKING) CODE
            % LFPmat(ch,:,trl) = ECoG2kp(ch,ceil(trialStart(trl)*Fs + win(1)*Fs):ceil(trialStart(trl)*Fs + win(2)*Fs));
            %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            %% constructing LFP tensor.
            LFPmat(:,tt,ch,aS) = dData(ch,round(trialStarts(tt)-pre)*params.Fs:round(trialStarts(tt)+post)*params.Fs);
            
            
            %% [20160610] calculating single trial spectrograms here.
            params.trialave = 0;
            %             if isequal(mod(tt,50),0)
            %                 display(sprintf('calculated spectrograms for for first %d trials out of %d...',tt,nTrials))
            %             end
            %             [S(:,:,tt),t,f] = mtspecgramc(LFPmat(:,tt,ch,aS),movingWin,params);
            
            
        end
        
        %         % calculating the median spectrograms for each conflict type.
        %         S0 = squeeze(nanmean(10*log10(S(:,:,trialType==1)),3));
        %         S1a = squeeze(nanmean(10*log10(S(:,:,trialType==2)),3));
        %         S1b = squeeze(nanmean(10*log10(S(:,:,trialType==3)),3));
        %         S2 = squeeze(nanmean(10*log10(S(:,:,trialType==4)),3));
        
        %         keyboard
        
        
        %% [20160610] calculating trial averaged spectrograms.
        params.trialave = 1;
        display('calculating trial averaged spectrogram for no conlfict trials...')
        [S0,~,~] = mtspecgramc(LFPmat(:,trialType==1,ch,aS),movingWin,params);
        display('calculating trial averaged spectrogram for simon conlfict trials...')
        [S1a,~,~] = mtspecgramc(LFPmat(:,trialType==2,ch,aS),movingWin,params);
        display('calculating trial averaged spectrogram for flanker conlfict trials...')
        [S1b,~,~] = mtspecgramc(LFPmat(:,trialType==3,ch,aS),movingWin,params);
        display('calculating trial averaged spectrogram for both conlfict trials...')
        [S2,t,f] = mtspecgramc(LFPmat(:,trialType==4,ch,aS),movingWin,params);
        
        
        %% time vectors
        tsec = linspace(-pre,post,size(LFPmat,1));
        tspec = linspace(-pre,post,length(t));
        
        % plotting LFP and spectrograms
        handleF = ch*100;
        brdr = 0.04;
        figure(handleF)
        % first, ERPs
        plotmultipleaxes(2,1,4,brdr,handleF)
        hold on
        plot(tsec,mean(LFPmat(:,trialType==1,ch,aS),2),'color',col0)
        plot(tsec,mean(LFPmat(:,trialType==2,ch,aS),2),'color',col1a)
        plot(tsec,mean(LFPmat(:,trialType==3,ch,aS),2),'color',col1b)
        plot(tsec,mean(LFPmat(:,trialType==4,ch,aS),2),'color',col2)
        hold off
        xlim([-1 2])
        set(gca,'linewidth',2,'fontsize',16)
        title(deblank(LFPlabels{ch}))
        
        % then, spectrograms.
        plotmultipleaxes(2,4,2,brdr,handleF)
        imagesc(tspec,f,normlogspec(S0)')
        xlim([-1 2])
        set(gca,'linewidth',2,'fontsize',16)
        axis xy square
        
        plotmultipleaxes(4,4,2,brdr,handleF)
        imagesc(tspec,f,normlogspec(S1a)')
        xlim([-1 2])
        set(gca,'linewidth',2,'fontsize',16)
        axis xy square
        
        plotmultipleaxes(6,4,2,brdr,handleF)
        imagesc(tspec,f,normlogspec(S1b)')
        xlim([-1 2])
        set(gca,'linewidth',2,'fontsize',16)
        axis xy square
        
        plotmultipleaxes(8,4,2,brdr,handleF)
        imagesc(tspec,f,normlogspec(S2)')
        xlim([-1 2])
        set(gca,'linewidth',2,'fontsize',16)
        axis xy square
        
        colorbar('NorthOutside')
        colormap(jet)
        
        
        %% saving figures.
            display('saving to file')
            dbPath = 'C:\Users\melete2\Desktop\TobyData\figs\LFP\';
            fName = sprintf('%s/%s_session_%d_Channel_%s_LFP_%sAligned',dbPath,patientID,sessionNum,deblank(LFPlabels{ch}),alignName);
            maximize(handleF)
            pause(3)
            saveas(handleF,fName,'pdf')
            close(handleF)
        
        
        display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
        
        %
        %             %% saving stats.
        %             if exist(['./' patientID],'dir')
        %                 try
        %                     fName = sprintf('./%s/Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
        %                     save([fName '.mat'],'neuronStats')
        %                 catch
        %                     mkdir(sprintf('./%s/Data/',patientID))
        %                     fName = sprintf('./%s/Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
        %                     save([fName '.mat'],'neuronStats')
        %                 end
        %             elseif exist('./Data','dir')
        %                 fName = sprintf('./Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
        %                 save([fName '.mat'],'neuronStats')
        %             else
        %                 fName = sprintf('%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
        %                 save([fName '.mat'],'neuronStats')
        %             end
        
    end % looping over lfp channels
end % looping over align spots (Stimulus & response)

