function dat = getNS5Data(dat,fn5,varargin)
p = inputParser;
p.addOptional('dsEye',30,@isnumeric);
p.addOptional('dsDiode',1,@isnumeric);
p.addOptional('nsEpoch',[0 0],@isnumeric);
p.parse(varargin{:});
downsampleeye = p.Results.dsEye;
downsamplediode = p.Results.dsDiode;
nsEpoch = p.Results.nsEpoch;


DIODE_CHAN = 3;
EYE_CHAN = [1 2 4]; % eye X, eye Y, pupil diameter
PUPIL_CHAN = 4;

hdr5 = read_nsx(fn5,'readdata',false);
ns5Samp = double(hdr5.hdr.Fs);
fprintf('Found %d channels of NS5 data.\n',hdr5.hdr.nChans);
%clockFs = double(data.hdr.clockFs);
for tind = 1:length(dat)
    epochStartTime = dat(tind).time(1) - nsEpoch(1);
    epochEndTime = dat(tind).time(2) + nsEpoch(2);
    nsEndTime = hdr5.hdr.nSamples / hdr5.hdr.Fs;
    if epochStartTime < 0
        epochStartTime = 0;
    end
    if epochEndTime > nsEndTime
        epochEndTime = nsEndTime;
    end
    msec = dat(tind).trialcodes(:,3);
    codes = dat(tind).trialcodes(:,2);
    codesamples = round(msec*ns5Samp);
    
    eyedata.codesamples = [codes codesamples];
    eyes = read_nsx(fn5,'chanindx',EYE_CHAN,'begsample',round(epochStartTime*ns5Samp),'endsample',round(epochEndTime*ns5Samp));
    eyedata.trial = downsample(eyes.data',downsampleeye)';
    eyedata.startsample = codesamples(1)/downsampleeye;
    eyedata.dataFs = ns5Samp/downsampleeye;
    
    diode.codesamples = [codes codesamples];
    diodes = read_nsx(fn5,'chanindx',DIODE_CHAN,'begsample',round(epochStartTime*ns5Samp),'endsample',round(epochEndTime*ns5Samp));
    diode.trial = int16(downsample(diodes.data,downsamplediode));
    diode.startsample = codesamples(1)/downsamplediode;
    diode.dataFs = ns5Samp/downsamplediode;
    
    pupil.codesamples = [codes codesamples];
    pupils = read_nsx(fn5,'chanindx',PUPIL_CHAN,'begsample',round(epochStartTime*ns5Samp),'endsample',round(epochEndTime*ns5Samp));
    pupil.trial = int16(downsample(pupils.data,downsamplediode));
    pupil.startsample = codesamples(1)/downsamplediode;
    pupil.dataFs = ns5Samp/downsamplediode;
    
    dat(tind).eyedata = eyedata;
    dat(tind).diode = diode;
    dat(tind).pupil = pupil;
    dat(tind).nsTime = (0:1:size(eyedata.trial,2)-1)./ns5Samp - nsEpoch(1);
end
end