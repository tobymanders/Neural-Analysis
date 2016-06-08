tic;
%load ('allUnits.mat');


%set constants
intensity=150;
tetrode=8;
ketamine=10080;
depth=29;
date='4-22-15';


s=1;
n = length(allUnits)+1;
allUnits(n).date=datenum(date);
allUnits(n).frequency=30000;
allUnits(n).region='ACC';
allUnits(n).animal='leonidas';
allUnits(n).spikes=A_001a;
allUnits(n).waveforms=A_001a_wf;
allUnits(n).trials=length(recFlick);
allUnits(n).latency=latency;
allUnits(n).depth=depth;
allUnits(n).laserOn=laserOn;
allUnits(n).recFlick=recFlick;
allUnits(n).PSTH=nex(:,s);
allUnits(n).intensity=intensity;
allUnits(n).tetrode=tetrode;
allUnits(n).ketamine=ketamine;
s=s+1;
% 
% n = length(allUnits)+1;
% allUnits(n).date=datenum(date);
% allUnits(n).frequency=30000;
% allUnits(n).region='ACC';
% allUnits(n).animal='leonidas';
% allUnits(n).spikes=A_001b;
% allUnits(n).waveforms=A_001b_wf;
% allUnits(n).trials=length(recFlick);
% allUnits(n).latency=latency;
% allUnits(n).depth=depth;
% allUnits(n).laserOn=laserOn;
% allUnits(n).recFlick=recFlick;
% allUnits(n).PSTH=nex(:,s);
% allUnits(n).intensity=intensity;
% allUnits(n).tetrode=tetrode;
% allUnits(n).ketamine=ketamine;
% s=s+1;
% 
% 
% n = length(allUnits)+1;
% allUnits(n).date=datenum(date);
% allUnits(n).frequency=30000;
% allUnits(n).region='ACC';
% allUnits(n).animal='leonidas';
% allUnits(n).spikes=A_001c;
% allUnits(n).waveforms=A_001c_wf;
% allUnits(n).trials=length(recFlick);
% allUnits(n).latency=latency;
% allUnits(n).depth=depth;
% allUnits(n).laserOn=laserOn;
% allUnits(n).recFlick=recFlick;
% allUnits(n).PSTH=nex(:,s);
% allUnits(n).intensity=intensity;
% allUnits(n).tetrode=tetrode;
% allUnits(n).ketamine=ketamine;
% s=s+1;

% 
% n = length(allUnits)+1;
% allUnits(n).date=datenum(date);
% allUnits(n).frequency=30000;
% allUnits(n).region='ACC';
% allUnits(n).animal='leonidas';
% allUnits(n).spikes=A_001d;
% allUnits(n).waveforms=A_001d_wf;
% allUnits(n).trials=length(recFlick);
% allUnits(n).latency=latency;
% allUnits(n).depth=depth;
% allUnits(n).laserOn=laserOn;
% allUnits(n).recFlick=recFlick;
% allUnits(n).PSTH=nex(:,s);
% allUnits(n).intensity=intensity;
% allUnits(n).tetrode=tetrode;
% allUnits(n).ketamine=ketamine;
% s=s+1;
% 
% 
% n = length(allUnits)+1;
% allUnits(n).date=datenum(date);
% allUnits(n).frequency=30000;
% allUnits(n).region='ACC';
% allUnits(n).animal='leonidas';
% allUnits(n).spikes=A_001e;
% allUnits(n).waveforms=A_001e_wf;
% allUnits(n).trials=length(recFlick);
% allUnits(n).latency=latency;
% allUnits(n).depth=depth;
% allUnits(n).laserOn=laserOn;
% allUnits(n).recFlick=recFlick;
% allUnits(n).PSTH=nex(:,s);
% allUnits(n).intensity=intensity;
% allUnits(n).tetrode=tetrode;
% allUnits(n).ketamine=ketamine;
% s=s+1;

%save ('allUnits.mat', 'allUnits','-v7.3');
% T=toc;
% load scriptTime;
% x=length(scriptTime);
% scriptTime(x+1,1)=n;
% scriptTime(x+1,2)=T;
% save('scriptTime.mat','scriptTime');