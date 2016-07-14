function [] = titleRAW(patientID)
% Appends 'Raw' to files in Raw folder
pathway = strcat('C:\Users\melete2\Desktop\TobyData\', num2str(patientID),'\Raw\');
names = dir(pathway);
names = names(3:end);
for y = 1:length(names)
    if ~strcmpi ('RAW', names(y).name(end-7:end-5))
        ext = names(y).name(end-3:end);
        movefile([pathway names(y).name], [pathway names(y).name(1:end-4) '_RAW_' ext]);
    end
end
