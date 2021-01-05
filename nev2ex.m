function [trials,ex] = nev2ex(fn, varargin)
%function ex = nev2ex3(fn)
%function ex = nev2ex3(fn, varargin)
%
% Function to read nev files created with Ex and convert to a
% structure with trial information.
%
% required argument:
% 
% fn: NEV file name, to read in 'nev', which contains
% a Nx3 matrix where the columns are elec num, sort code, time(s)
%
% optional arguments:
%
% alignCode: if present, adjusts the times so that a time of 0
%  corresponds to the presence of the sort code given.  If more than one 
%  instance of alignCode is found in the trial, the first is used.
%
% collapseConditions: if true, will throw all conditions together into one
%  cell (i.e., units X 1 X repeats). Default is false.
%
% keepTrialCode: defaults to 5 (REWARD), the output will only include
%  trials that have this code. If a '1' is input then all trials are returned.
%
% readLFP - reads 1kHz LFP data from NS2 file
% readEyes - reads Eye X, Y and Pupil Diameter from NS5 file
% readDiode - reads 30 kHz photodiode signal from NS5 file
%  * All of these should be logicals
%  * Checks if NS2/NS5 file exist and errors if they don't
%  * Assumes channels 1-4 contain eye and diode data, and downsamples them
%  from 30kHz to 1kHz.
%
% convertEyes - if true, converts eye X/Y values to degrees
%
% nsEpoch - 2-element vector, with amount of time in seconds to pad each
% trial with. If [1 2] is passed, then each trial's NS data will have an
% extra 1 s of samples before the '1' code and 2 s of samples after the
% '255' code.
%
% output is a struct with the fields:
%   EVENTS: spike times, in seconds, unit x condition x sort code
%   CODES: digital codes and timestamps (s), unit x condition x
%   sort code
%   CHANNELS: List of units, elec num x sort code
%   AND MORE ...
%

% Could add a "verbose" mode?

%
% OLD STUFF
% Diode alignment now disabled!!
% WARNING: DIODE FEATURE NOT WELL TESTED!!!!!!!
% diode: the 30kHz diode signal (optional). should call read_nsx and then
% send in the diode signal (usually channel 3) to this function
%
% diodeThresh: the value for the diode threshold.  aligning happens the
%  first time the diode is greater than this value. If thresh is not
%  specified, it uses 75% of the max value in the diode signal


%% optional input arguments
p = inputParser;
p.addOptional('readDiode',false,@islogical);
p.addOptional('readLFP',false,@islogical);
p.addOptional('readEyes',false,@islogical);
p.addOptional('convertEyes',false,@islogical);
p.addOptional('alignCode',1,@isscalar);
p.addOptional('keepTrialCode',[],@isscalar);
p.addOptional('collapseConditions',false,@islogical);
p.addOptional('nsEpoch',[0 0],@isnumeric);

p.parse(varargin{:});

readDiode = p.Results.readDiode;
readLFP = p.Results.readLFP;
readEyes = p.Results.readEyes;
convertEyes = p.Results.convertEyes;
alignCode = p.Results.alignCode;
keepTrialCode = p.Results.keepTrialCode;
collapseConditions = p.Results.collapseConditions;
nsEpoch = p.Results.nsEpoch;


%% save some important values
START_TRIAL = 1;
END_TRIAL = 255;

% warn if multiple instances of the align code are found in a trial
warnAlignFlag = 0;

% flag to include 0 and 255 codes in output. By default, remove these
include_0_255 = 0;

% number of trials to read before printing message
readerStatusInterval = 100;

% Channel defaults
DIODE_CHAN = 3;
EYE_CHAN = [1 2 4]; % eye X, eye Y, pupil diameter

% setup a default keep code. Could use '1' here if the default is to keep
% all trials
if isempty(keepTrialCode)
    keepTrialCode = 5; % this is the REWARD code, could also use CORRECT which is 150
end

if alignCode ~= 1
    alignFlag = 1;
else
    alignFlag = 0;
end

%% Read in the NEV data and NSX if requested
nev = readNEV(fn);
[pathstr,fnt]=fileparts(fn);
fn2 = [pathstr,'/',fnt,'.ns2'];
fn5 = [pathstr,'/',fnt,'.ns5'];
nsEndTime = NaN;

if readLFP
    if exist(fn2,'file')~=2
        error(['Unable to find file - ',fn2]);
    end
    hdr2 = read_nsx(fn2,'readdata',false);
    lfpSamp = hdr2.hdr.Fs;
    fprintf('Found %d channels of LFP data.\n',hdr2.hdr.nChans);
    if size(hdr2.hdr.timeStamps,2) > 1
        error('File was paused - not supported');
    end
    lfpChan = hdr2.hdr.label;
    lfpChan = str2double(lfpChan);
    nsEndTime = hdr2.hdr.nSamples / hdr2.hdr.Fs;
end

if readDiode || readEyes
    if exist(fn5,'file')~=2
        error(['Unable to find file - ',fn5]);
    end
    hdr5 = read_nsx(fn5,'readdata',false);
    ns5Samp = hdr5.hdr.Fs;
    fprintf('Found %d channels of NS5 data.\n',hdr5.hdr.nChans);
    if hdr5.hdr.nChans < (numel(EYE_CHAN) + numel(DIODE_CHAN))
        error('Not enough channels found in NS5 file');
    end
    if ~strcmp(hdr5.hdr.label{1},'10241')
        % channels 10241-10244 are the usual ones, anything else 
        % means some other data may be in this NS5
        error('First channel must be 10241 - there are probably unusual data streams in this NS5 file');
    end
    nsEndTime = hdr5.hdr.nSamples / hdr5.hdr.Fs;
end

%% start parsing the NEV data

codes = nev(nev(:,1)==0,2:3);
nev = nev(nev(:,1) ~= 0,:);

channels = unique(nev(:,1:2),'rows');
if ~include_0_255
  channels = channels(channels(:,2) ~= 0 & channels(:,2) ~= 255 & channels(:,1) ~= 0,:);
else
  channels = channels(channels(:,1) ~= 0,:);
end

starts = find(codes(:,1) == START_TRIAL);
ends = find(codes(:,1) == END_TRIAL);

% warn if different number of trial starts and ends, but we deal with this
% below to handle interrupted trials
if length(starts) ~= length(ends)
    warning([num2str(length(starts)),' Trial starts and ',num2str(length(ends)),' Trial ends. Continuing.']);
end

%% fill up a cell array of data from each trial
trials = cell(length(starts),1);

tstarts = starts;
tcount = 1;
fprintf(' Completed trials: \n');
while ~isempty(tstarts) % look through all the start codes

    if rem(tcount,readerStatusInterval)==0
        fprintf(' %i ... \n',tcount);
    end
    
    % find the first end code after the start you're looking at
    nextend = ends(find((ends - tstarts(1)) > 0, 1 ));
    tcodes = codes(tstarts(1):nextend,:);

    % Double-check there's no more than one end code
    if length(find(tcodes(:,1) == 255)) > 1
        error(['Trial ',num2str(tcount),' had more than 1 end code. Serious error - exiting']);
    end
    
    % check to see if there is an interrupted trial buried here -
    % if so, cut off the codes after the second '1' in the trial so
    % that you can parse just this trial
    if length(find(tcodes(:,1) == 1)) ~= 1
        warning(['Trial ',num2str(tcount),' was apparently interrupted.']);
        % if this is the last trial, exit the loop
        if (length(tstarts)==1)
            trials(end)=[];
            tstarts(1) = []; % delete this trial from the list
            continue
        else % if it's got a buried 1, trim off from there on and keep this interrupted trial
            tcodes = tcodes(1:find(tcodes(:,1)==1,1,'last')-1,:);
            nextend = tstarts(1)+size(tcodes,1);
        end
    end
    
    trial = struct();
    trial.codes = tcodes;
    cndIndex = find(trial.codes(:,1)>32768,1);

    if collapseConditions % throw all conditions into one cell, or separate them
        trial.cnd = 1;
    else
        trial.cnd = trial.codes(cndIndex,1)-32768;
    end
    trial.codes(cndIndex,:) = [];
    
    startTime = trial.codes(1,2);
    endTime = trial.codes(end,2);    
    trial.startTime = startTime;
    trial.endTime = endTime;

    % MATT 03/30/17 - I don't understand how this error check could ever fail
    % Check codes by start/endtime and compare to indexing method
    tc=codes(codes(:,2) >= startTime & codes(:,2) <= endTime,:);
    if (size(tc,1) ~= size(tcodes,1))
        error(['Code timing vs. indexing mismatch on trial ',num2str(tcount)]);
    end
    
    trial.spikes = nev(nev(:,3) >= startTime & nev(:,3) < endTime,:);

    % adjust times to start of trial
    trial.spikes(:,3) = trial.spikes(:,3) - startTime;
    trial.codes(:,2) = trial.codes(:,2) - startTime;

    % extra buffer on beginning and end of NSX data, but make sure it stays
    % within the bounds of the NSX file
    epochStartTime = startTime - nsEpoch(1);
    epochEndTime = endTime + nsEpoch(2);
    
    if epochStartTime < 0
        epochStartTime = 0;
    end
    if epochEndTime > nsEndTime
        epochEndTime = nsEndTime;
    end

    if readDiode
        diode = read_nsx(fn5,'chanindx',DIODE_CHAN,'begsample',round(epochStartTime*ns5Samp),'endsample',round(epochEndTime*ns5Samp));
        trial.diode = diode.data(1:ns5Samp/lfpSamp:end); %downsample
        trial.nsTime = (0:1:size(trial.diode,2)-1)./double(lfpSamp)  - nsEpoch(1);
    end
    
    if readEyes
        eyes = read_nsx(fn5,'chanindx',EYE_CHAN,'begsample',round(epochStartTime*ns5Samp),'endsample',round(epochEndTime*ns5Samp));
        trial.eyes = eyes.data(:,1:ns5Samp/lfpSamp:end); %downsample
        trial.nsTime = (0:1:size(trial.eyes,2)-1)./double(lfpSamp)  - nsEpoch(1);
    end
    
    if readLFP
        lfp = read_nsx(fn2,'begsample',round(epochStartTime*lfpSamp),'endsample',round(epochEndTime*lfpSamp));
        trial.lfp = lfp.data;
        trial.nsTime = (0:1:size(trial.lfp,2)-1)./double(lfpSamp) - nsEpoch(1); 
    end

    trial.channels = channels;
    msgInd = trial.codes(:,1) >= 256 & trial.codes(:,1) < 512;
    trial.msgs = char(trial.codes(msgInd,1)-256)';
    variables = regexp(trial.msgs,';','split');
    for j = 1:length(variables)
        if variables{j}
            try
                % this option is for plain numeric values
                eval(['trial.env.' variables{j} ';']);
            catch
                % if it's not simple numeric, save it as a string
                vt=regexp(variables{j},'=','split');
                try
                    eval(['trial.env.' char(vt{1}(:))' '=' '''' char(vt{2}(:))' '''' ';']);
                    %eval(['trial.env.' char(vt{1}(1)) '=' '''' char(vt{1}(2)) '''' ';']);
                catch
                    % if that fails, send a warning 
                    %if ~warnVarFlag
                    disp(['*** Could not parse variable name: ',variables{j}]);
                    %disp(['No more warnings of this type will be displayed']);
                    %   warnVarFlag = 1;
                    %end
                end
            end
        end
    end
    
    trials{tcount} = trial;
    tcount = tcount + 1;
    tstarts(1) = []; % delete this trial from the list
end
fprintf('\n');

%% remove invalid trials
trials = cell2mat(trials);
cndlist = arrayfun(@(x) x.cnd,trials);
cnds = unique(cndlist);

% find trials with the keepTrialCode (usually '5' for REWARD)
goodtrials = arrayfun(@(x) sum(x.codes(:,1) == keepTrialCode)>0,trials);

% find trials with both a start and an end - this could be optional if we 
% wanted to rescue them somehow. For now we just delete the incomplete trials
started = arrayfun(@(x) sum(x.codes(:,1) == START_TRIAL)>0,trials);
ended = arrayfun(@(x) sum(x.codes(:,1) == END_TRIAL)>0,trials);
completed = started + ended - 1;
completed(completed<0) = 0;

% only include completed trials in the good trials list
goodtrials = boolean(goodtrials .* completed);

%% Fill up cell arrays of trial data
EVENTS = cell(size(channels,1),max(cnds),1);
MSGS = cell(max(cnds),1);
NSTIME = cell(max(cnds),1);
if readDiode
    DIODE = cell(max(cnds),1);
end
if readEyes
    EYES = cell(numel(EYE_CHAN),max(cnds),1);
end
if readLFP
    LFP = cell(numel(lfpChan),max(cnds),1);
end
CODES = cell(max(cnds),1);
ENV = cell(max(cnds),1);

% find all the codes sent between trials and put them in 'params'
% This is wildly inefficient, it's a kludge to make sure you can reference
% the parameters easily from each trial. We could improve this - MATT
preTrial = cell(0);
lastEnd = 0;
params = struct();
for i = 1:length(trials)
    tri = trials(i);
    preCodes = codes(codes(:,2) > lastEnd & codes(:,2) < tri.startTime,1);
    preTrial{i} = [char(preCodes(preCodes >= 256 & preCodes < 512) - 256)'];
    if ~isempty(preCodes)
%         disp(['Trial # ',num2str(i),' preceded by digital codes']);
        %disp(preTrial{i});
    end
    variables = regexp(preTrial{i},';','split');
    for j = 1:length(variables)
        k = strfind(variables{j},'=');
        if k            
            lhs = variables{j}(1:k-1);
            rhs = variables{j}(k+1:end);
            % MATT - should we check here to see if any value has
            % changed and report back if it has?
            try
                eval(['params.' lhs '=[' rhs '];']);
            catch
                eval(['params.' lhs '=''' rhs ''';']);
            end
        end
    end
    % put all the params into the trials struct
    trials(i).env = catstruct(trials(i).env, params);    
    
    lastEnd = tri.endTime;
end

% This is a stupid hack to make eye2deg work later
params.block = params;

%% Align each trial's spikes and codes, then store in the cell arrays
for i = 1:max(cnds)
    
    % status message
    disp(['Converting trials for condition ',num2str(i),' of ',num2str(max(cnds))]);
    
    theseTrials = trials(cndlist==i & goodtrials);
    
    for j = 1:length(theseTrials)
        if alignFlag
            if alignCode < 0 % alignCode -1 was used to align to diode
                error('NEV2EX no longer supports aligning on the diode');
            else
                if (numel(alignCode)==1)
                    alignTime = theseTrials(j).codes(theseTrials(j).codes(:,1) == alignCode,2);
                elseif (numel(alignCode)>=1)
                    codestr=num2str(theseTrials(j).codes(:,1)');
                    patstr = [];
                    for I=1:length(alignCode)-1
                        patstr = [patstr,num2str(alignCode(I)),' \s '];
                    end
                    patstr = [patstr,num2str(alignCode(I+1))];
                    idx = regexp(codestr,patstr,'start');
                    alignTime = 0;
                end
                if isempty(alignTime)
                    fprintf(['Repeat %i of condition %i does not ' ...
                                  'have align code %i'],j,i, ...
                                 alignCode);
                    alignTime = 0;
                elseif length(alignTime) > 1
                    if ~warnAlignFlag
                        fprintf(['Repeat %i of condition %i has ' ...
                                      '%i occurrences of align code %i ' ...
                                      '- using 1st occurrence'],j,i, ...
                                     length(alignTime),alignCode);
                        disp('No more warnings of this type will be displayed');
                        warnAlignFlag = 1;
                    end
                    alignTime = alignTime(1);
                end
            end
        else
            alignTime = 0;
        end % of alignment

        % put the spikes into the right spot
        for k = 1:size(channels,1)
            valid = theseTrials(j).spikes(:,1) == channels(k,1) ...
                        & theseTrials(j).spikes(:,2) == channels(k,2);
                 
            EVENTS{k,i,j} = theseTrials(j).spikes(valid,3) - alignTime;
        end

        CODES{i,j} = theseTrials(j).codes;
        CODES{i,j}(:,2) = CODES{i,j}(:,2) - alignTime;
        MSGS{i,j} = theseTrials(j).msgs;
        ENV{i,j} = theseTrials(j).env;
        
        if readDiode
            DIODE{i,j} = theseTrials(j).diode;
        end
        if readEyes
            eyedat = theseTrials(j).eyes;
            if convertEyes % make the X/Y in deg vis angle
                eyedat(1:2,:) = eye2deg(eyedat(1:2,:),params);
            end
            for k = 1:numel(EYE_CHAN)
                EYES{k,i,j} = eyedat(k,:);
            end
        end
        if readLFP
            for k = 1:numel(lfpChan)
                LFP{k,i,j} = theseTrials(j).lfp(k,:);
            end
        end
        
        if readDiode || readEyes || readLFP
            NSTIME{i,j} = theseTrials(j).nsTime-alignTime;            
        end
    end
end

%% Put all the data into the ex struct for output
ex = struct();
ex.EVENTS = EVENTS;
if readDiode
    ex.DIODE = DIODE;
end
if readEyes
    ex.EYES = EYES;
end
if readLFP
    ex.LFP = LFP;
    ex.LFPCHANNELS = lfpChan;
end
if readDiode || readEyes || readLFP
    ex.NSTIME = NSTIME;
end

ex.MSGS = MSGS;
ex.CODES = CODES;
ex.CHANNELS = channels;
ex.TRIAL_SEQUENCE = cndlist(goodtrials);
ex.REPEATS = hist(ex.TRIAL_SEQUENCE,1:max(cnds));
ex.PRE_TRIAL = preTrial(goodtrials);
ex.ENV = ENV;
