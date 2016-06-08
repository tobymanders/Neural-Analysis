% load RHD2000 file
tic;
read_Intan_RHD2000_file;

% get sampling rate
sample_rate = frequency_parameters.amplifier_sample_rate;

% determine number of channels
n = length(amplifier_channels);

% create and sort table of tetrode values
channelMap = xlsread('C:\Users\Toby\Documents\Toby Data\channelMappingforIntanandAdapter');
trodeNumber = zeros(n,2);
for y=1:n
    trodeNumber(y,1)=channelMap(find(channelMap(:,3)==amplifier_channels(:,y).native_order),2);
    trodeNumber(y,2)=amplifier_channels(:,y).native_order;
end

% bandpass filter amplifier data
Wn = [300/(0.5*sample_rate) 7500/(0.5*sample_rate)];

[b,a] = butter (2, Wn);
for m = 1:n
    amplifier_data(m,:) = filtfilt(b,a,amplifier_data(m,:));
end

% write tetrodes
ext = '.nex';
tet = 'tetrode';
for l=1:8
    trodeChan = find (trodeNumber(:,1)>=l & trodeNumber(:,1)<(l+1));
    if ~isempty(trodeChan)
        nexFile = nexCreateFileData (sample_rate);
        a = length(trodeChan);
        for x=1:a
            nexFile = nexAddContinuous (nexFile, 0, sample_rate, amplifier_data(trodeChan(x),:), amplifier_channels(:,trodeChan(x)).native_channel_name);
        end
        writeNexFile (nexFile, strcat (path, tet, num2str(l), ext));
    end
    clearvars trodeChan
end

% write laser events to excel file
difference = diff(board_dig_in_data);
test = find (difference==1);
time = test/sample_rate;
time = transpose(time);
xlswrite (strcat (path, 'laserOn.xlsx'), time);

T=toc;
clearvars -except T filesize;
load('generateFilesTime.mat');
generateFilesTime(end+1,1)=filesize;
generateFilesTime(end,2)=T;
save('generateFilesTime.mat','generateFilesTime');
fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', T);
clear