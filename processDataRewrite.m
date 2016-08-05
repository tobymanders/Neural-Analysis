function [] = processData(patientID, region, varargin)

% Filters data for spike and LFP analysis, extracts hippocampal channels,
% detects (and optionally removes) discharges and noise, and saves data as
% .nex and .mat for further analysis/sorting.
%
% First input must be patient ID number, e.g. '11'. Second argument is
% either can be 'hippocampus,' 'acc' or 'all'.
%
% Optional arguments:
%
%   'LFPonly':              Does not perform tasks specific to spike analysis.
%
%   'detectDischarges':     Detects discharges and saves variable in both .mat
%                           file and .nex file for sorting/later analysis
%
%   'removeDischarges':     After detecting discharges, will remove cues
%                           and spikes from discharge periods. By default, cues will be removed
%                           within 5 seconds of beginning or end of discharge. Note: this option
%                           will enable detectDischarge, too.
%
%   'plot':                 Plots results of denoising process. Off by
%                           default.
%
%   'spartan':              Only performs filters for LFP and spikes and
%                           saves results.
tic;
memSample(1) = getfield(memory,'MemUsedMATLAB');
%Rename files to append '_RAW_'
titleRAW(patientID);

%% default parameters
lfponly = 0; % Disables spike steps when enabled.
dischargeFlag = 0; %Enables discharge detection
removeDischarge = 0; %Strips discharge cues from LFP data
preDisch = 5;
postDisch = 5; %Specifies amount of time before and after discharges in which to eliminate cues
filepath = 'C:\Users\ShethLab\Desktop\TobyData\';
cd ([filepath num2str(patientID) '\Raw\'])
files = dir ('*.ns5');
files = {files.name}';
plotFlag = 0;

%% parse inputs
for i = 1:length(varargin)
    inputArgument = varargin{i};
    if strcmpi(inputArgument, 'LFPonly')
        lfponly = 1;
    elseif strcmpi (inputArgument, 'detectDischarges')
        dischargeFlag = 1;
    elseif strcmpi (inputArgument, 'removeDischarges')
        dischargeFlag = 1;
        removeDischarge = 1;
    elseif strcmpi (inputArgument, 'plot')
        plotFlag = 1;
    elseif strcmpi (inputArgument, 'spartan')
        spartanFlag = 1;
        plotFlag = 0;
        dischargeFlag = 0;
        removeDischarge = 0;
        lfponly = 0;
    end
end


%% Loop over all .ns5 files in Raw directory
for curr = 1:length(files)
    file = char(files(curr));
    %% Open NSx file
    openNSx (strcat(filepath, num2str(patientID), '\Raw\', file),'skipfactor',10);
    
    %% Find channels with hippocampal data and find recording frequency
    namesAll = {NS5.ElectrodesInfo.Label};
    if strcmpi (region, 'hippocampus')
        regCh = find ((strncmp ('uHH',namesAll,3) == 1) | (strncmp ('uHT',namesAll,3) == 1)...
        | (strncmp ('uHB',namesAll,3) == 1) | (strncmp ('uHC',namesAll,3) == 1)| (strncmp ('uLHC',namesAll,3) == 1)...
        | (strncmp ('uRHC',namesAll,3) == 1) | (strncmp ('HH',namesAll,2) == 1) |(strncmp ('HT',namesAll,2) == 1)...
        |(strncmp ('HC',namesAll,2) == 1)|(strncmp ('LHC',namesAll,3) == 1));
    elseif strcmpi (region, 'acc')
        regCh = find ((strncmp ('uAC',namesAll,3) == 1));
    elseif strcmpi (region, 'all')
        regCh = find ((strncmp ('ainp',namesAll,4) ~= 1));
    end
    sample_rate = NS5.MetaTags.SamplingFreq;
    samples = 10*length(NS5.Data);
    memSample(end+1) = getfield(memory,'MemUsedMATLAB');
    deci = double(NS5.Data(regCh,:));
    names = namesAll(regCh);
    channelNum = length(regCh);
    clear NS5;

    
    if ~isempty (regCh)
        %% Find discharges
        if dischargeFlag == 1
            deciAvg = mean(deci);
            avg = mean(deciAvg);
            stdev = std(deciAvg);
            plot (abs(deciAvg));
            hold;
            discharges = find(abs(deciAvg) > (2*stdev + avg));
        end
        
        %% Load channels with data
        options = sprintf('c:%d:%d',min(regCh),max(regCh));
        NS5 = openNSx(strcat(filepath, num2str(patientID), '\Raw\', file), options, 'p:double');
        
        %% Find noise
        %decideci = deciAvg(1:10:end);
        raw = nan(channelNum,samples);
        raw = NS5.Data;
        clear NS5;
        memSample(end+1) = getfield(memory,'MemUsedMATLAB');
        
        %% Bandpass filter amplifier data for spikes
        Wn = [300/(0.5*sample_rate) 7500/(0.5*sample_rate)];
        [b,a] = butter (4, Wn);
        highPass = nan(size(raw));
        for ch = 1:channelNum
            highPass(ch,:) = filtfilt(b,a,raw(ch,:));
        end    
            
        %deciData = nan(channelNum,length(raw)/10); %move this line
        save(strcat(filepath,num2str(patientID),file(1:end-8),upper(region),'_PRESORT'), 'raw', 'highPass', '-v7.3');
        memSample(end+1) = getfield(memory,'MemUsedMATLAB');
        
        %% Low pass filter amplifier data for LFP
        Wn = 300/(0.5*sample_rate);
        [b,a] = butter (2, Wn);
        lowPass = nan(size(raw));
        for ch = 1:channelNum
            lowPass(ch,:) = filtfilt(b,a,raw(ch,:));
        end 
        save(strcat(filepath,num2str(patientID),file(1:end-8),upper(region),'_PRESORT'),'-append', 'lowPass');
        memSample(end+1) = getfield(memory,'MemUsedMATLAB');
        clear lowPass
        
        %% Bandpass filter amplifier data for discharge detection
        if dischargeFlag == 1
            Wn = [5/(0.5*sample_rate) 7500/(0.5*sample_rate)];
            [b,a] = butter (2, Wn);
            disPass = filtfilt(b,a,raw(1,:)')';
            
            %% Detect discharge periods
            smoothCross = abs(smooth(hippoStruct(1).highFive,10000));
            dischMean = mean(smoothCross);
            dischStd = std(smoothCross);
            cutoff = dischMean + 2*dischStd;
            idx1 = smoothCross>cutoff;
            idx1(1) = 0;
            idx = find(idx1);
            yest = smoothCross(idx-1)<cutoff;
            posCross = idx(yest);
            posCross(end) = [];
            idx(end) = 0;
            yest2 = smoothCross(idx+1)<cutoff;
            negCross = idx(yest2);
            discharges = zeros (length(posCross),2);
            for a = 1:length(posCross)
                discharges(a,1) = posCross(a);
                later = find(negCross>posCross(a));
                discharges(a,2) = negCross(later(1));
            end
            memSample(end+1) = getfield(memory,'MemUsedMATLAB');
            
        end
        %% Find and remove noisy epochs
        if spartanFlag ~= 1
            noise = mean(deciData);
            noisePower = power(noise,2);
            smoothPower = smooth(noisePower, 30000);
            SPA = mean(smoothPower);
            SPstd = std(smoothPower);
            cutoff = SPA + SPstd;
            locations = 10*find(smoothPower>cutoff);
            filledIn = zeros(10*length(locations),1);
            for o = 1:(length(filledIn)-1)
                filledIn(o) = locations(floor(o/10)+1)+(rem(o,10)-1);
            end
            filledIn(end)=[];
            memSample(end+1) = getfield(memory,'MemUsedMATLAB');
            
            %remove noisy epochs from all channels
            for q = 1:channelNum
                for p = 1:length(filledIn)
                    hippoStruct(q).denoisedData(filledIn(p))=0;
                    LFP(q,filledIn(p))=0;
                end
            end
            
            %% Results
            if plotFlag == 1;
                if lfponly ~=1
                    channel = 1; %select channel to view results from
                    
                    plot0 = subplot(4,1,1);
                    plot(hippoStruct(channel).data);
                    title('Raw Data');
                    
                    filtPlot = subplot(4,1,2);
                    plot(hippoStruct(channel).filteredData);
                    title('Bandpassed Data');
                    
                    plot1 = subplot(4,1,3);
                    plot (smoothPower, 'LineWidth', 2);
                    title('Smoothed Data Denoising');
                    
                    plot2 = subplot (4,1,4);
                    plot (hippoStruct(channel).denoisedData);
                    title('Final Denoised Data');
                    deciNoise = filledIn(1:1000:length(filledIn));
                    for num = 1:length(deciNoise);
                        SP = deciNoise(num);
                        line([SP SP],get(plot2,'YLim'),'Color',[1 0 0], 'LineWidth', 0.1, 'LineStyle', '--');
                    end
                    
                    print(strcat(filepath, num2str(patientID),'\Presort\', 'Results_', '_', upper(region), '_', file(1:end-9),'_PRESORT_'),'-dpng');
                end
            end
        end
        %% Get Elliot's behavioral data
        NEV = openNEV(strcat(filepath, num2str(patientID), '\Raw\', file(1:end-4),'.nev'), 'nomat', 'nosave');
        trigs = double(NEV.Data.SerialDigitalIO.UnparsedData);
        trigTimes = double(NEV.Data.SerialDigitalIO.TimeStampSec);
        TimeRes = NEV.MetaTags.TimeRes;
        memSample(end+1) = getfield(memory,'MemUsedMATLAB');

        nTrials = sum(trigs>199 & trigs <207);
        %% parsing behavior
        trialType = zeros(1,nTrials);
        sub = ones(nTrials,1);
        condition = trigs(find(trigs>199 & trigs<207) - 2*sub);
        
        
        %% These are the correct codes. Double Checked on 20160216
        trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
        trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
        trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
        trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)
        
        
        %% Cues and response times.
        
        cueTimes = trigTimes(find(trigs>199 & trigs<207) - 2*sub);
        respTimes = trigTimes(find(trigs>199 & trigs<207) - sub);
        
        if length(cueTimes) == length(respTimes);
            RTs = respTimes - cueTimes;
        end
        
%         if (length(respTimes) < length(cueTimes))
%             cueTimes(end) = [];
%         end
%         
%         if length(trialType)>length(respTimes)
%             trialType(end) = [];
%         end
%         
%         cueTimes(find ((respTimes-cueTimes)>20,1,'first')) = [];
%         
%         if length(respTimes)>length(cueTimes)
%             respTimes(end) = [];
%         end
%         
        if spartanFlag ~= 1
            for b = 1:length(cueTimes)
                if find(cueTimes(b)==filledIn) > 0
                    cueTimes(b) = 0;
                    respTimes(b) = 0;
                    if removeDischarge == 1
                        if any((cueTimes(b)+postDisch)>discharges(:,2) & cueTimes(b-preDisch)<discharges(:,1))
                            cueTimes(b) = 0;
                            respTimes(b) = 0;
                        end
                    end
                end
            end
        end
        
        
        
%         RTs = respTimes-cueTimes;
        memSample(end+1) = getfield(memory,'MemUsedMATLAB');

        
        %% Locate easy, intermediate and hard trials %new
        easy = find (trialType == 1);
        cueEasy = cueTimes(easy);
        respEasy = respTimes(easy);
        
        int1 = find (trialType == 2);
        cueInt1 = cueTimes(int1);
        respInt1 = respTimes(int1);
        
        int2 = find (trialType == 3);
        cueInt2 = cueTimes(int2);
        respInt2 = respTimes(int2);
        
        hard = find (trialType == 4);
        cueHard = cueTimes(hard);
        respHard = respTimes(hard);
        
        % Remove cues from noisy sections
        
        cueTimes(cueTimes == 0) = [];
        respTimes(respTimes == 0) = [];
        RTs (RTs == 0) = [];
        cueEasy(cueEasy == 0) = [];
        cueInt1(cueInt1 == 0) = [];
        cueInt2(cueInt2 == 0) = [];
        cueHard(cueHard == 0) = [];
        
        
        %% do statistics on RTs over conflict
        [P,~,~] = anova1(RTs,trialType,'off');
        
        %% plot RTs over conflict
        figure
        boxplot(RTs,trialType)
        title([ ', omnibus p-value = ' num2str(P)])
        set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'none','spatial','distractor','both'})
        xlabel('Conflict Type')
        ylabel('RT (seconds)')
        ylim([0 5])
        axis square
        
        
        %% write data to nex file
        destination = strcat(filepath, num2str(patientID),'\Presort\');
        
        ext = '.nex';
        nexFile = nexCreateFileData (sample_rate);
        clearvars NS5
        for l=1:channelNum
            if lfponly ~= 1
                chName = names{l};
                nexFile = nexAddContinuous (nexFile, 0, sample_rate, highPass(l,:), chName(1:4));
            end
            nexFile = nexAddEvent (nexFile, cueTimes, 'cueTimes');
            nexFile = nexAddEvent (nexFile, respTimes, 'respTimes');
            nexFile = nexAddEvent (nexFile, cueEasy, 'cueEasy');
            nexFile = nexAddEvent (nexFile, cueInt1, 'cueInt1');
            nexFile = nexAddEvent (nexFile, cueInt2, 'cueInt2');
            nexFile = nexAddEvent (nexFile, cueHard, 'cueHard');
            nexFile = nexAddEvent (nexFile, respEasy, 'respEasy');
            nexFile = nexAddEvent (nexFile, respInt1, 'respInt1');
            nexFile = nexAddEvent (nexFile, respInt2, 'respInt2');
            nexFile = nexAddEvent (nexFile, respHard, 'respHard');
            if dischargeFlag ~= 0
                nexFile = nexAddEvent (nexFile, posCross, 'discStart');
                nexFile = nexAddEvent (nexFile, negCross, 'discStop');
            end
        end
        if lfponly ~= 1
            writeNexFile (nexFile, strcat (destination, file(1:end-8),'_PRESORT_',upper(region),'_', ext));
        else if lfponly == 1
                writeNexFile (nexFile, strcat (destination, file(1:end-8),'_LFP_',upper(region),'_', ext));
            end
        end
        clearvars nexFile
        ext = '.mat';
        if lfponly ~= 1 && dischargeFlag == 1
            save (strcat(destination, file(1:end-8),'_PRESORT_',upper(region),'_', ext), 'discharges','LFP', '-v7.3')
        else if lfponly == 1 && dischargeFlag == 1
                save (strcat(destination, file(1:end-8),'_LFP_',upper(region),'_', ext), 'LFP', 'discharges', '-v7.3')
                else if lfponly == 1 && dischargeFlag ~= 1
                        save (strcat(destination, file(1:end-8),'_LFP_',upper(region),'_', ext), 'LFP', '-v7.3')
                    end
                end
            end
        end
   %If there is no data from that region
    end
    clearvars  LFP discharges deciData
%Loop over files
% Print time & plot memory usage
plot(memSample);
title('Memory Usage');
clc;
% clearvars -except T;
T = toc;
fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', T);