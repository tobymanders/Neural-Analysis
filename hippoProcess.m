%% Finds hippocampal data, extracts and filters it, removes noisy sections, and saves as a .nex file for further analysis.

%% Open NSx file
tic;
openNSxNew
file = NS5.MetaTags.Filename(1:end-4);
filepath = NS5.MetaTags.FilePath(1:end-3);

%% Find channels with hippocampal data and find recording frequency
names = {NS5.ElectrodesInfo.Label};
hippoCh = find (strncmp ('uHH',names,3) == 1);
sample_rate = NS5.MetaTags.SamplingFreq;
channelNum = length(hippoCh);


%% Save hippocampal data to a new structure
hippoStruct(channelNum).name=[];
hippoStruct(channelNum).channel=[];
hippoStruct(channelNum).data=[];
for x = 1:channelNum
    hippoStruct(x).name = names(hippoCh(x));
    hippoStruct(x).channel = hippoCh(x);
    hippoStruct(x).data = double(NS5.Data(hippoCh(x),:));
end


%% Bandpass filter amplifier data for spikes
Wn = [300/(0.5*sample_rate) 7500/(0.5*sample_rate)];

[b,a] = butter (2, Wn);
for r = 1:channelNum
    hippoStruct(r).filteredData = filtfilt(b,a,hippoStruct(r).data);
    deciData(r,:) = hippoStruct(r).filteredData(1:10:end);
    hippoStruct(r).denoisedData = hippoStruct(r).filteredData;
    
end

%% Low pass filter amplifier data for LFP
Wn = 300/(0.5*sample_rate);

[b,a] = butter (2, Wn);
for r = 1:channelNum
    hippoStruct(r).LFP = filtfilt(b,a,hippoStruct(r).data);
end

%% Bandpass filter amplifier data for discharge detection
Wn = [5/(0.5*sample_rate) 7500/(0.5*sample_rate)];

[b,a] = butter (2, Wn);
for r = 1:channelNum
    hippoStruct(r).highFive = filtfilt(b,a,hippoStruct(r).data);    
end

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
yest2 = smoothCross(idx+1)<cutoff;
negCross = idx(yest2);
discharges = zeros (length(posCross),2);
for a = 1:length(posCross)
    discharges(a,1) = posCross(a);
    later = find(negCross>posCross(a));
    discharges(a,2) = negCross(later(1));
end
hippoStruct(1).dischargePeriods = discharges;

%% Find and remove noisy epochs
comb = zeros (1,length(deciData));
for n = 1:length(deciData(1,:))
    comb(n) = mean (deciData(:,n));
end
noisePower = (power((deciData(1,:)-deciData(2,:)-deciData(3,:)-deciData(4,:)-deciData(5,:)-deciData(6,:)-deciData(7,:)-deciData(8,:)),2));
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

%remove noisy epochs from all channels
for q = 1:channelNum
    for p = 1:length(filledIn)
        hippoStruct(q).denoisedData(filledIn(p))=0;
    end
end

%% Results

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

print(strcat(filepath,'Presort\','Results_',file),'-dpng');

%% Get Elliot's behavioral data
NEV = openNEV(strcat(filepath,'Raw\', file,'.nev'));
trigs = double(NEV.Data.SerialDigitalIO.UnparsedData);
trigTimes = double(NEV.Data.SerialDigitalIO.TimeStampSec);
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<28);


%% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% Cues and response times.
cueTimes = trigTimes(trigs>=1 & trigs<=27);
respTimes = trigTimes(trigs>=100 & trigs<=103);
if (length(cueTimes) > length(respTimes) && length(cueTimes) == length(trialType))
    cueTimes(end) = [];
    trialType(end) = [];
end

RTs = respTimes-cueTimes;

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

%% do statistics on RTs over conflict
[P,~,behavioralStats] = anova1(RTs,trialType,'off');

%% plot RTs over conflict
sessionNum = 13;
figure
boxplot(RTs,trialType)
title([ ', omnibus p-value = ' num2str(P)])
set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'none','spatial','distractor','both'})
xlabel('Conflict Type')
ylabel('RT (seconds)')
ylim([0 5])
axis square


%% write data to nex file
ext = '.nex';
destination = strcat(filepath,'Presort\');
nexFile = nexCreateFileData (sample_rate);
clearvars NS5
for l=1:channelNum
    chName = names{hippoCh(l)};
    nexFile = nexAddContinuous (nexFile, 0, sample_rate, hippoStruct(l).denoisedData, chName(1:4));
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
    nexFile = nexAddEvent (nexFile, posCross, 'discStart');
    nexFile = nexAddEvent (nexFile, negCross, 'discStop');


end
writeNexFile (nexFile, strcat (destination, file, ext));
clearvars nexFile
ext = '.mat';
save (strcat(destination, file, ext), 'hippoStruct', 'discharges', '-v7.3')

% Print time
clc;
clearvars -except T;
T = toc;
fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', T);
