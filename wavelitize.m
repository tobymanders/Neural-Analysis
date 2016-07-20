function [ POWER, PHASE_CONSISTENCY] = wavelitize( srate, numcycles, num_freqs, T1, T2, B1, B2, dats, dims)

%T1 = index(time start)
%T2 = index(time end)
%B1 = index(baseline time start)
%B2 = index(baseline time end)
%dats = chan x time x data
%dims = size(dats)
%srate = sampling rate
%numcycles = number of cycles (ex: 4)
%numfreqs = max freq


%srate=EEG.srate; numcycles=4; num_freqs=50;
%     if strcmp(condition,'cues')
%         B1=find(tx==-2000);
%         B2=find(tx==0);
%     elseif strcmp(condition,'resp')
%         B1=find(tx==-2000);
%         B2=find(tx==-500);
%     end
% 
% T1=find(tx==-2000);
% T2=length(tx);


frex=logspace(.01,1.7,num_freqs); s=numcycles./(2*pi.*frex); t=-2:1/srate:2;

    for fi=1:length(frex)
        w(fi,:) = exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*s(fi)^2));
    end

DAT4STUFF=reshape(squeeze(dats(chan,:,:)),1,dims(2)*dims(3)); %concatenate epochs

    for fi=1:num_freqs
        stuffnjunk=fconv_JFC(DAT4STUFF,w(fi,:));  % use fconv_JFC for fast convolution
        stuffnjunk=stuffnjunk((size(w,2)-1)/2+1:end-(size(w,2)-1)/2);
        stuffnjunk=reshape(stuffnjunk,dims(2),dims(3));

        BASE = abs(stuffnjunk(B1:B2,:)).^2;
        %

            POWER(fi,:) = mean(abs(stuffnjunk(T1:T2,:)).^2,2);
            % dB Conversion
            POWER(fi,:) = 10*(log10(squeeze(POWER(fi,:))) - log10(repmat(mean(BASE),1,size(tx2disp,2))) ); 
            PHASE_CONSISTENCY(fi,:) = abs(mean(exp(1i*(  angle(stuffnjunk(T1:T2,:))  )),2));  % Derives to Phase Locking Value (Lachaux et al.)


        clear stuffnjunk BASE;
    end

end

