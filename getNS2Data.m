function dat = getNS2Data(dat,fn2,varargin)
p = inputParser;
p.addOptional('nsEpoch',[0 0],@isnumeric);
p.parse(varargin{:});
nsEpoch = p.Results.nsEpoch;

hdr2 = read_nsx(fn2,'readdata',false);
ns2Samp = double(hdr2.hdr.Fs);
fprintf('Found %d channels of LFP data.\n',hdr2.hdr.nChans);
%clockFs = double(data.hdr.clockFs);
for tind = 1:length(dat)
    epochStartTime = dat(tind).time(1) - nsEpoch(1);
    epochEndTime = dat(tind).time(2) + nsEpoch(2);
    nsEndTime = hdr2.hdr.nSamples / hdr2.hdr.Fs;
    if epochStartTime < 0
        epochStartTime = 0;
    end
    if epochEndTime > nsEndTime
        epochEndTime = nsEndTime;
    end
    msec = dat(tind).trialcodes(:,3);
    codes = dat(tind).trialcodes(:,2);
    codesamples = round(msec*ns2Samp);
    
    lfpdata.codesamples = [codes codesamples];
    lfp = read_nsx(fn2,'begsample',round(epochStartTime*ns2Samp),'endsample',round(epochEndTime*ns2Samp));
    lfpdata.trial = lfp.data;
    lfpdata.startsample = codesamples(1);
    lfpdata.dataFs = ns2Samp;
    lfpChan = hdr2.hdr.label;
    lfpdata.chan = str2double(lfpChan);
    dat(tind).lfp = lfpdata;
    dat(tind).nsTime = (0:1:size(lfpdata.trial,2)-1)./ns2Samp - nsEpoch(1);
end
end