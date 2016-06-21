function [] = collectMSITnevs(patientID)
%findMSITnevs This function scours the patient folder structure to find
%MSIT data.
%   Input should be a string with the patient ID, e.g. 'CUBF11'.
filepath = '/mnt/mfs/patients/BF/';
cd([filepath patientID]);
files = subdir ('*.nev');
dest = [filepath 'MSITnevs'];
for x = 1:length(files)
    NEV = openNEV ([files(x).name],'nosave', 'noread', 'nowarning');
    trigs = double(NEV.Data.SerialDigitalIO.UnparsedData);
    if sum(trigs==90) > 100
        copyfile (files(x).name, dest);
    end
end