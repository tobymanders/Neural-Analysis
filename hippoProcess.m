%% Finds hippocampal data, extracts and filters it, removes noisy sections, and saves as a .nex file for further analysis.

%% Open NSx file
tic;
openNSxNew

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


%% Bandpass filter amplifier data
Wn = [300/(0.5*sample_rate) 7500/(0.5*sample_rate)];

[b,a] = butter (2, Wn);
for m = 1:channelNum
    hippoStruct(m).filteredData = filtfilt(b,a,hippoStruct(m).data);
end

%% Find and remove noise sections from data

for m = 1:channelNum
    hippoStruct(m).denoised = hippoStruct(m).filteredData;
    deciData = hippoStruct(m).filteredData(1:10:end);
    for x = 2:length(deciData)
        absDiff(x) = abs(deciData(x)-deciData(x-1));
    end
    smoothDiff = abs(smooth(absDiff, 30000));
    meanVal = median(smoothDiff);
    s = std(smoothDiff);
    index = find(smoothDiff>(meanVal+2*s));
    for x = 2:length(index)
        index(x,2) = index(x) - index(x-1);
    end
    borders = find (index(:,2) > 30);
    for x = 1:(length(borders)-1)
        if median(smoothDiff(index((borders(x):borders(x+1)),1)))>meanVal
            for l = borders(x):borders(x+1)
                hippoStruct(m).denoised(l*10) = 0;
            end
        else
        end
    end
    noiseGauge(1,length(smoothDiff)*10) = 0;
    for y = 1:(length(borders)-1);
        start   = index(borders(y),1);
        stop    = index(borders(y+1)-1,1);
        sectionLength = (stop-start);
        subData = smoothDiff(start:stop);
        if median(subData)>(meanVal+s)
            for p = (10*start):(10*stop)
                noiseGauge(p)=1;
                hippoStruct(m).denoised(p)=0;
            end
        end
    end
    
    for z = 2:length(smoothDiff);
        noiseEdge(z)=(noiseGauge(z)~=noiseGauge(z-1));
    end
    edges = find(noiseEdge==1);
    
    %% Results
   
    plot0 = subplot(4,1,1);
    plot(hippoStruct(m).data);
    title('Raw Data');
    
    filtPlot = subplot(4,1,2);
    plot(hippoStruct(m).filteredData);
    title('Bandpassed Data');
    
    plot1 = subplot(4,1,3);
    plot (smoothDiff, 'LineWidth', 2);
    for num = 1:length(edges);
        SP = edges(num);
        line([SP SP],get(plot1,'YLim'),'Color',[1 0 0], 'LineWidth', 0.25, 'LineStyle', '--');
    end
    title('Smoothed Data Denoising');
    
    hline1 = refline (0, meanVal+2*s);
    hline1.Color = 'r';
    hline1.LineStyle = '--';
    hline2 = refline (0, meanVal-2*s);
    hline2.Color = 'r';
    hline2.LineStyle = '--';
    subplot (4,1,4);
    plot (hippoStruct(m).denoised);
    title('Final Denoised Data');
end




% write data to nex file
ext = '.nex';
destination = 'C:\Users\melete2\Desktop\TobyNex\';
nexFile = nexCreateFileData (sample_rate);
file = NS5.MetaTags.Filename(1:end-4);
for l=1:channelNum
    chName = names{hippoCh(l)};
    nexFile = nexAddContinuous (nexFile, 0, sample_rate, hippoStruct(l).denoised, chName(1:4));
end
writeNexFile (nexFile, strcat (destination, file, ext));

% Print time
T = toc;
fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', T);