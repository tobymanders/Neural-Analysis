function [] = collectMSITnevs(patientID)
%Finds MSIT data and copies .nev and .ns5 files to new folder
%   Searches all subdirectories in patient directory for recordings with
%   >100 trials of MSIT data, and copies relevant data for analysis.
filepath = '/mnt/mfs/patients/BF/';
cd([filepath patientID]);
files = subdir ('*.nev');
dest = '/mnt/mfs/home/NIMASTER/ehs2149/TobyData/MSITnevs';
n = 0;
for x = 1:length(files)
    if files(x).bytes >40000
        NEV = openNEV ([files(x).name],'nosave', 'noread', 'nowarning');
        trigs = double(NEV.Data.SerialDigitalIO.UnparsedData);
        if sum(trigs==90) > 100
            copyfile (files(x).name, dest);
            copyfile ([files(x).name(1:end-3) 'ns5'], dest);
            n = n+1;
        end
    else
    end
end
display(sprintf ('Found and moved %d recordings', n));
