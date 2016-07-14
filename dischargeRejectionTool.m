function [rejectTrials] = dischargeRejectionTool(patientID,sessionNum,nevFile,plotFlag)
%DISCHARGREJECTIONTOOL Tool for rejecting trials with epileptiform discharges.
%
%   [rejectTrials] = dischargeRejectionTool(patientID,sessionNum,nevFile)
%   rejects outliers based on their maximum and mimimum values.
%
%   optional input argument: plotFlag == 0 will not plot rejected trials.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)

% default argument for plotting.
if nargin==3
    plotFlag = 1;
end


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
Fs = 2e3;


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


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


%% NS3.Data => dData (rather than de-noising)
dData = double(NS3.Data);


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% initializing LFPmat
LFPlabels = {NS3.ElectrodesInfo.Label};
LFPmat = zeros(((pre+post)*Fs)+1,nTrials,length(LFPlabels),2);
% for loop to save multiple epochs
for aS = 1
    %% going to start by just aligning on CUE.
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
        
        if plotFlag
            % figure window.
            ah = figure('Color',[0 0 0]);
            hold on
        end
        
        for tt = 1:nTrials
            
            
            %% [20160622] constructing LFP tensor.
            LFPmat(:,tt,ch,aS) = dData(ch,round(trialStarts(tt)-pre)*Fs:round(trialStarts(tt)+post)*Fs);
            
            
            %% [20160622] time vectors
            tsec = linspace(-pre,post,size(LFPmat,1));
            
            if plotFlag
                %% [20160622] plotting LFP from each trial brushing to reject
                plot(tsec,LFPmat(:,tt,ch,aS)+((tt-1)*1000),'color',rgb('white'))
                axis tight off
            end
            
        end
        rejectTrials(ch).rejectTheseTrials = outliers(squeeze(range(LFPmat(:,:,ch,aS))));
        rejectTrials(ch).channelLabel = deblank(LFPlabels{ch});
        
        %%
        if plotFlag
            for rjct = rejectTrials(ch).rejectTheseTrials
                %% [20160622] plotting LFP from each trial brushing to reject
                plot(tsec,LFPmat(:,rjct,ch,aS)+((rjct-1)*1000),'color',rgb('springgreen'))
            end
            maximize(gcf)
            title(sprintf('Channel: %s. || Rejected Trials in Green',rejectTrials(ch).channelLabel),'color',rgb('springgreen'))
            
            
            %% saving figures.
            try
                display('saving to Elliot"s dropbox. Thanks!')
                dbPath = sprintf('/home/elliot/Dropbox/MSITunits_emu/%s/Figs',patientID);
                fName = sprintf('%s/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrials(ch).channelLabel);
                saveas(gcf,fName, 'pdf')
                close(gcf)
            catch
                maximize(ch*1000+un)
                if exist(['../../' patientID],'dir')
                    try
                        fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrials(ch).channelLabel);
                        saveas(gcf,fName, 'pdf')
                        close(gcf)
                    catch
                        mkdir(sprintf('../../%s/Figs/',patientID))
                        fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrials(ch).channelLabel);
                        saveas(gcf,fName, 'pdf')
                        close(gcf)
                    end
                elseif exist('../../Figs','dir')
                    fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrials(ch).channelLabel);
                    saveas(gcf,fName, 'pdf')
                    close(gcf)
                else
                    fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrials(ch).channelLabel);
                    saveas(gcf,fName, 'pdf')
                    close(gcf)
                end
            end
            
            display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
        end
    end % looping over lfp channels
    hold off
end % looping over align spots (Stimulus & response)

