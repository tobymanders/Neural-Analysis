% Extract hippocampal data, bandpass filter, and save as .nex

% Open NSx file
tic;
openNSxNew

% Find channels with hippocampal data and find recording frequency
names = {NS5.ElectrodesInfo.Label};
hippoCh = find (strncmp ('uHH',names,3) == 1);
sample_rate = NS5.MetaTags.SamplingFreq;
channelNum = length(hippoCh);

% Save hippocampal data to a new structure and delete other variables
hippoStruct(channelNum).name=[];
hippoStruct(channelNum).channel=[];
hippoStruct(channelNum).data=[];
for x = 1:channelNum
    hippoStruct(x).name = names(hippoCh(x));
    hippoStruct(x).channel = hippoCh(x);
    hippoStruct(x).data = double(NS5.Data(hippoCh(x),:));
end

% bandpass filter amplifier data
Wn = [300/(0.5*sample_rate) 7500/(0.5*sample_rate)];

[b,a] = butter (2, Wn);
for m = 1:channelNum
    hippoStruct(m).filteredData = filtfilt(b,a,hippoStruct(m).data);
end


% write data to nex file
ext = '.nex';
destination = 'C:\Users\melete2\Desktop\TobyNex\';
nexFile = nexCreateFileData (sample_rate);
file = NS5.MetaTags.Filename(1:end-4);
for l=1:channelNum
    chName = names{hippoCh(l)};
    nexFile = nexAddContinuous (nexFile, 0, sample_rate, hippoStruct(l).filteredData, chName(1:4));
end
writeNexFile (nexFile, strcat (destination, file, ext));

% Print time
T = toc;
fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', T);