%% plotPSTH

function y = plotPSTH(patientID,sessionNum)

%% Pull data
NEX = readNexFile;
events = [NEX.events{1:20,1}];
units = [NEX.neurons{:}];
names =  {units.name};
%
% unsorted = find(~cellfun('isempty',regexpi (names, 'U$')'));
% unita = find(~cellfun('isempty',regexpi (names, 'a$')'));
% unitb = find(~cellfun('isempty',regexpi (names, 'b$')'));
% unitc = find(~cellfun('isempty',regexpi (names, 'c$')'));

%% Create neural timing variable
ChanUnitTimestamp = [];
for m = 1:length(units)
    chID = units(m).wireNumber + 1;
    uID = units(m).unitNumber +1;
    dim = length(units(m).timestamps);
    ChanUnitTimestamp = [ChanUnitTimestamp; (ones(dim,1)*chID) ones(dim,1)*uID units(m).timestamps];
end


%% Pull data from nev file
[nevf, nevp] = uigetfile;
[~, nvName, ~] = fileparts([nevp nevf]);
NEV = openNEV([nevp nevf], 'nomat', 'nosave');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs>199 & trigs<207);
% detect and remove false starts or practice
if max(diff(trigTimes) > 6)
    display(sprintf('\nIt seems there were practice trials in this file...\n\nRemoving %d events before event time gap.',find(diff(trigTimes) > 6)))
    trigs(1:find(diff(trigTimes) > 6)) = [];
    trigTimes(1:find(diff(trigTimes) > 6)) = [];
    nTrials = sum(trigs>199 & trigs<207);
end


%% timing (seconds)
pre = 2;
post = 3;


%% creating neural timing variable

inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);
nTrials = length(events(1).timestamps);


%% aligning AP data.
display('Aligning AP data on stimulus and response.');
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
    
    
    % looping over Channels
    for ch = 1:nChans
        
        % looping over number of units in the AP data
        nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
        for un = 1:nUnits
            
            % getting unit times for the current channel and unit.
            unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
            
            % loooping over trials
            for tt = 1:nTrials
                
                %% putting the data in a structure
                data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post) - repmat(trialStarts(tt)-pre,length(unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
                
            end
        end
    end
    
    % save data structure?
    saveFlag = 0;
    if saveFlag
        display(['saving spike time structure for ' alignName '...'])
        save([patientID '_session' num2str(sessionNum) '_spikeTimeStruct_alignedon' alignName '.mat'],'data','pre','post','alignName')
    end
    
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



%% Rasters and PSTHs
RT = events(2).timestamps - events(1).timestamps;

% Alternate ranking of raster by response time
for rankRT = 1:2
    switch rankRT
        case 1
            [~,I] = sort(RT);
            rankFlag = 'RT';
        case 2
            I = [1:nTrials]';
            rankFlag = 'Chronology';
    end
    
    for aS2 = 1:length(data)
        % which alignment spot
        switch aS2
            case 1
                alignName = 'Cue';
                trialStarts =  events(1).timestamps;
            case 2
                alignName = 'Response';
                trialStarts =  events(2).timestamps;
        end
        
        % looping over units
        for ch = 1:size(data(aS2).channel,2)
            for un = 2:size(data(aS2).channel(ch).unit,2)
                
                display('plotting raster over conflict... grab a cocktail, this may take a while.')
                %% plotting rasters and PSTHs
                figure(aS2)
                ah_ras = plotmultipleaxes(1,1,2,0.08,aS2);
                hold on
                for tt = 1:nTrials
                    % changing raster color based on trial type
                    if trialType(I(tt))==1
                        rasCol = rgb('DarkGreen');
                    elseif trialType(I(tt))==2
                        rasCol = rgb('Goldenrod');
                    elseif trialType(I(tt))==3
                        rasCol = rgb('OrangeRed');
                    elseif trialType(I(tt))==4
                        rasCol = rgb('FireBrick');
                    end
                    
                    % plotting rasters for conflict (in the least efficient way possible)
               
                    for sp = 1:size(data(aS2).channel(ch).unit(un).trial(I(tt)).times,1)
                        try
                            line([data(aS2).channel(ch).unit(un).trial(I(tt)).times(sp)-pre data(aS2).channel(ch).unit(un).trial(I(tt)).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', rasCol)
                        catch
                            line([data(aS2).channel(ch).unit(un).trial(I(tt)).times(sp)-pre data(aS2).channel(ch).unit(un).trial(I(tt)).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', 'k')
                        end
                    end
                    
                    % Add response/cue time line
                    switch aS2
                        case 1
                            line([RT(I(tt)) RT(I(tt))], [tt-(9/20) tt+(9/20)],'linewidth',4, 'color', 'Black')
                        case 2
                            line([-RT(I(tt)) -RT(I(tt))], [tt-(9/20) tt+(9/20)],'linewidth',4, 'color', 'Black')
                            
                    end

                    
                    % stimulus timing lines
                    line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
                    % raster plot details
                    xlim([-1 2])
                    ylim([0 nTrials])
                    str = sprintf('patient %s, Channel %d, Unit %d; aligned on %s; sorted by %s',patientID ,ch ,un-1,alignName, rankFlag);
                    title(str,'fontsize',18);
                    ylabel('Trials','fontsize', 16)
                    set(gca, 'linewidth', 2, 'fontsize', 16);
                                        
                    
                    
                end
                hold off
                
                
                %% PSTH=
                % calculating psths
                kernelWidth = 25  ./1000;
                [Reasy,t,Eeasy] = psth(data(aS2).channel(ch).unit(un).trial(trialType(1:nTrials)==1), kernelWidth, 'n', [0 pre+post]);
                [Rhard,t,Ehard] = psth(data(aS2).channel(ch).unit(un).trial(trialType(1:nTrials)==4), kernelWidth, 'n', [0 pre+post]);
                
                
                try
                    [Rmedi1,t,Emedi1] = psth(data(aS2).channel(ch).unit(un).trial(trialType(1:nTrials)==2), kernelWidth, 'n', [0 pre+post]);
                catch
                    display('no spatial conflict')
                end
                
                try
                    [Rmedi2,t,Emedi2] = psth(data(aS2).channel(ch).unit(un).trial(trialType(1:nTrials)==3), kernelWidth, 'n', [0 pre+post]);
                catch
                    display('no distracter conflict')
                end
                
                
                tsec = t-repmat(pre,1,length(t));
                
                
                %% generate statistics for each Neuron
                [pEasy,Heasy] = ranksum(Reasy(tsec>=-1 & tsec<=-0.5),Reasy(tsec>=0 & tsec<=1.5));
                [pHard,Hhard] = ranksum(Rhard(tsec>=-1 & tsec<=-0.5),Rhard(tsec>=0 & tsec<=1.5));
                if isequal(aS,1)
                    neuronStats.Cue{ch,un} = table(pEasy,pHard,'VariableNames',{'pEasy','pHard'});
                elseif isequal(aS,2)
                    neuronStats.Response{ch,un} = table(pEasy,pHard,'VariableNames',{'pEasy','pHard'});
                end
                
                
                %% make it pretty
                % colors and colors and colors and colors
                EcolEasy = rgb('darkgreen');
                EcolHard = rgb('darkred');
                EcolMedi1 = rgb('Orangered');
                EcolMedi2 = rgb('goldenrod');
                
                % PSTH plot
                ah_psth = plotmultipleaxes(2,1,2,0.08,aS2);
                hold on
                
                
                %% plotting psths
                try
                    % plotting med psth
                    patch([tsec fliplr(tsec)],[Rmedi1+Emedi1 fliplr(Rmedi1-Emedi1)], EcolMedi1,'edgecolor','none','facealpha',0.5)
                    plot(tsec,Rmedi1,'color',rgb('orangered'),'linewidth',2)
                catch
                    display('no spatial conflict')
                end
                
                try
                    % plotting med psth
                    patch([tsec fliplr(tsec)],[Rmedi2+Emedi2 fliplr(Rmedi2-Emedi2)], EcolMedi2,'edgecolor','none','facealpha',0.5)
                    plot(tsec,Rmedi2,'color',rgb('goldenrod'),'linewidth',2)
                catch
                    display('no spatial conflict')
                end
                
                % plotting easy psth
                patch([tsec fliplr(tsec)],[Reasy+Eeasy fliplr(Reasy-Eeasy)], EcolEasy,'edgecolor','none','facealpha',0.5)
                plot(tsec,Reasy,'color',rgb('DarkGreen'),'linewidth',2)
                
                % plotting hard psth
                patch([tsec fliplr(tsec)],[Rhard+Ehard fliplr(Rhard-Ehard)], EcolHard,'edgecolor','none','facealpha',0.5)
                plot(tsec,Rhard,'color',rgb('darkRed'),'linewidth',2)
                
                % stimulus timing lines
                line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
                
                % PSTH plot details
                xlim([-1 2])
                ylim([0 max(Rhard+Ehard)+3])
                hold off
                xlabel('Time (seconds)', 'fontsize', 16);
                ylabel('Firing Rate (spikes/second)', 'fontsize', 16);
                set(gca, 'linewidth', 2, 'fontsize', 16);
                
                
                %% saving figures.
                fpath = 'C:\Users\ShethLab\Desktop\TobyData\figures\';
                ss = (sprintf('%s_session_%d_Channel_%d_Unit_%d_Conflict_%s_aligned_%s_ordered',patientID,sessionNum,inclChans(ch),un-1,alignName,rankFlag));
                fName = [fpath ss];
                saveas(aS2,fName, 'pdf')
                close (aS2);
                
                %% saving stats.
                th = 'C:\Users\melete2\Desktop\TobyData\';
                ss = (sprintf('%s_session_%d_Channel_%d_Unit_%d_Conflict_%s_aligned_%s_ordered',patientID,sessionNum,inclChans(ch),un-1,alignName, rankFlag));
                fName = [fpath ss];
                save([fName '.mat'],'neuronStats')
                
            end
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot.
end % looping over align spots (Stimulus & response)


