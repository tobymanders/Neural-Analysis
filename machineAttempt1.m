NEX = readNexFile;
prepost = .3; %in seconds. Looks before and after.
tally = 0;
events = [NEX.events{1:10,1}];
neurons = [NEX.neurons{1:end,1}];
for k = 1:4
    eventtimes = events(k+6).timestamps; %easy, int1, int2, hard (response times)
    for j = 1:length(neurons)
        stamps = neurons(j).timestamps;
        baseline = length(stamps)/max(stamps);
        for i = 1:length(eventtimes)
            begin = eventtimes(i) - prepost;
            close = eventtimes(i) + prepost;
            unitfreq(i+tally,j) = (sum(stamps > begin & stamps < close) / (2*prepost))/baseline; %relative frequency compared with baseline
            trialtype{i+tally,1} = events(k+6).name;
        end
    end
    tally = length(unitfreq);
    clear eventtimes stamps baseline
end
A = array2table(unitfreq);
B = table(categorical(trialtype));
C = [A B];