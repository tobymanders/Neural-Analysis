%% Get file locations
[nexf, nexp] = uigetfile;
matf = [nexf(1:end-4) '.mat'];
MATOBJ = matfile([nexp matf]);
[~,nChan] = size(MATOBJ, 'hippoStruct');
% nChan = 8;
LFP = MATOBJ.LFP;
LFP = LFP(:,1:30:end);

%% Pull events from .nex file
NEX = readNexFile([nexp nexf]);
events = [NEX.events{1:10,1}];
nTrials = length(events(1).timestamps);
% 
% for y = 1:nChan
%     LFPlabels{y,1} = NEX.contvars{y,1}.name;
% end

clearvars NEX;


%% Pull LFP data from .mat file
% display('Extracting LFP data...')
% 
% for h = 1:nChan
%     tempStruct = MATOBJ.hippoStruct(1,h);
%     LFP(h,:) = tempStruct.LFP(1:30:end);
%     clearvars tempStruct
%     display(sprintf ('%d percent complete...',round((h-1)*100/nChan)))
% end


%% defining Chronux parameters.
movingWin = [0.8 0.080];
params.Fs = 1e3; % sampling frequency for LFP
params.pad = 2; %
params.fpass = [0 200]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 0; % average over trials {CHANGES BELOW}
params.err = [2 0.01]; % population error bars


%% timing (seconds)
pre = 2;
post = 3;


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


%% building LFP tensor and plotting ERPs and spectrograms.


LFPmat = zeros(((pre+post)*params.Fs)+1,nTrials,nChan,2);


% for loop to save multiple epochs
for aS = 1:2
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  events(1).timestamps;
        case 2
            alignName = 'Response';
            trialStarts =  events(2).timestamps;
    end
    
    for ch = 1:nChan
        for tt = 1:nTrials
            
            
            
            %% constructing LFP tensor.
            LFPmat(:,tt,ch,aS) = LFP(ch,round(trialStarts(tt)-pre)*params.Fs:round(trialStarts(tt)+post)*params.Fs);
            
            
            %% [20160610] calculating single trial spectrograms here.
            params.trialave = 0;
            
        end
        
        
%% parsing behavior
trialType = zeros(1,nTrials);
[~,easyCues] = ismember (events(3).timestamps, events(1).timestamps);
[~,int1Cues] = ismember (events(4).timestamps, events(1).timestamps);
[~,int2Cues] = ismember (events(5).timestamps, events(1).timestamps);
[~,hardCues] = ismember (events(6).timestamps, events(1).timestamps);


%% setting up codes for PSTHs over conflict types.
% These are the correct codes. Double Checked on 20160216
trialType(easyCues) = 1;    % Type 0 (Cond # 1-3)
trialType(hardCues) = 4;   % Type 2 (Cond # 4-15)
trialType(int1Cues) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(int2Cues) = 3;  % Type 1b Distractor interference (Cond # 21-27)

        
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
%         title(deblank(LFPlabels{ch}))
        
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
        
        
        
    end % looping over lfp channels
end % looping over align spots (Stimulus & response)

