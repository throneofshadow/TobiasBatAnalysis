function findcalls_session_tobias(wd,fs,varargin)
%% Find, extract, and save social vocalizations from audio recordings
% INPUT:
% wd: working directory, file path which contains files to be analyzed
% fs: sampling rate of recorded files
% Optional: 
% fileType: 'wav' or 'mat' to indicate format of recordings
% anal_dir: file path to save extracted files
%
% WARNING: Many individual parameters to be tuned!
%
% OUTPUT:
% Individual calls will be saved off in individual file in anal_dir named 
% with the original files's name and appeneded with a string '_Call_XXX'
% where XXX is the three digit number of calls within that file.

shift = 1; % is this post-shift? 1 for yes my friend.
%Collect unique identifiers from wd
[date_array,subj_array,time_array]=unique_gen(wd);
date_array = date_array';
subj_array = subj_array';
time_array = time_array';
if isempty(date_array )
    date_array = dir([wd, '*.mat']);
    subj_array = 0;
    time_array=0;
end
% Now I run findcalls over each unique date,time,subj combo inside the wd
for ggx = 1:length(date_array)
    Date = date_array(ggx);
    ExpTime = time_array(ggx);
    Subj = subj_array(ggx);
    Date =char(Date);
    ExpTime = char(ExpTime);
    Subj = char(Subj);
    % all of these parameters are now inside the wd loop!
    pnames = {'fileType', 'outputDir', 'dataVarName','findcall_inputs'};
    dflts  = {'wav', fullfile(wd, 'Analyzed_auto'),'recsGroup',{}};
    [fileType,outputDir,dataVarName,findcall_inputs] = internal.stats.parseArgs(pnames,dflts,varargin{:});
    % modified here to only pull in rec_files with Date,Subj,ExpTime
    rec_files = dir([wd,'*',Subj,'*',Date,'*',ExpTime,['*.' fileType]]);

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    n_files = length(rec_files);
    % pull in parameters and trigger info
    pfileType = '.txt.';
    paramTypes = '_param';
    param_file = dir([wd, '*',Date,'*',ExpTime,'*', paramTypes,'*',pfileType]);
    naming_rights = string(param_file.name);
    if isempty(param_file)
        % manually set params here
        lowpass = 0;
        highpass = 0;
        sessID = 0;
        sessionType = 0;
        
    elseif contains(naming_rights,'RecOnly')
        paramOpen = fopen(fullfile(wd,param_file.name));
        paramOpen = textscan(paramOpen,'%s','delimiter','\n');
        paramOpen = paramOpen{1};
        paramOpen = strjoin(paramOpen);
        paramOpen = string(paramOpen);
        lowpass = extractBetween(paramOpen,"high-pass filter frequency threshold:","Hz");
        highpass = extractBetween(paramOpen,"low-pass filter frequency threshold:", "Hz");
        sessID = extractBetween(paramOpen,'Session #:', '.');
        if isempty(sessID)
            sessID = 0;
        end
        sessionType = extractBetween(paramOpen,'Session type:','.');
        if isempty(sessionType)
            sessionType = 'cpAD';
        end
        lowcheck = contains(lowpass,'None'); % designed to check 
        highcheck = contains(highpass,'None');
        if lowcheck
            lowpass = 0;
        end
        if highcheck
            highpass = 0;
        end
    else
        paramOpen = fopen(fullfile(wd,param_file.name));
        paramOpen = textscan(paramOpen,'%s','delimiter','\n');
        paramOpen = paramOpen{1};
        paramOpen = strjoin(paramOpen);
        paramOpen = string(paramOpen);
        lowpass = extractBetween(paramOpen,"high-pass filter frequency threshold:","Hz");
        highpass = extractBetween(paramOpen,"low-pass filter frequency threshold:", "Hz");
        sessID = extractBetween(paramOpen,'Session #:', '.');
        if isempty(sessID)
            sessID = 0;
        end
        sessionType = extractBetween(paramOpen,'Session type:','.');
        if isempty(sessionType)
            sessionType = 'cpAD';
        end
        lowcheck = contains(lowpass,'None'); % designed to check 
        highcheck = contains(highpass,'None');
        if lowcheck
            lowpass = 0;
        end
        if highcheck
            highpass = 0;
        end
     
    
    %%% parse in triggers
    DataFileStruc = dir(fullfile(wd, sprintf('%s_%s_%s*events.txt', Subj, Date, ExpTime)));
    if isempty(DataFileStruc)
        FullStamps = [];
    end
    if contains((DataFileStruc.name),'RecOnly')
        FullStamps = [];
    else
        Fid_Data = fopen(fullfile(DataFileStruc.folder,DataFileStruc.name));
        EventsHeader = textscan(Fid_Data, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n',1);
        for hh=1:length(EventsHeader)
            if strfind(EventsHeader{hh}{1}, 'SampleStamp')
                EventsStampCol = hh;
            elseif strfind(EventsHeader{hh}{1}, 'Type')
                EventsEventTypeCol = hh;
            end
        end
        Events = textscan(Fid_Data, '%s\t%f\t%s\t%s\t%f\t%f\t%f');
        fclose(Fid_Data);
        VocId = find(strcmp('Vocalization', Events{EventsEventTypeCol}));
        FullStamps = Events{EventsStampCol}(VocId);
        %handle for no params or voc_triggers
        % Now I have stamps for the triggered times.
        %%%Generate hyperparams for data from param file

        % paramgen returns params from the param file
        % paramgen also returns the trigger timestamps
        % triggerstamps saved as FullStamps
    end
    for fln = 1:n_files
        filename=rec_files(fln).name;
        disp(filename)
        disp(['Analyzing file: ' num2str(fln) ' of ' num2str(n_files)])
        disp('...')
        % cut and load in parameters to file
        % edited for Tobias
        % parse filename for info && assign to hyperparams
        split_filename = strsplit(filename,'_'); 
        callNum = 1; % parameters to adjust
        batName = split_filename{1};
        saveDate = split_filename{2};
        saveTime = split_filename{3};
        callType = split_filename{4};
        micType = split_filename{5};
        recNum=extractBefore(split_filename{6},'.wav');
        % load data and calculate envelope
        data_raw = load_audio_data(wd,filename,fileType,dataVarName); 
        % Now, I need to filter using ba filter
        [b,a,zed] = butter(3,400/(fs/2),'high');
        [sos,g] = zp2sos(b,a,zed);
        senv = filtfilt(sos,g,data_raw);
        % Filter End
        filename = filename(1:end-4);
        if fln < n_files
            next_filename = rec_files(fln+1).name;
        else
            next_filename = [];
        end
        findcall_inputs = {findcall_inputs{:} 'filename' filename 'next_filename' next_filename};
        % Add in call for struct hyper_params to send in to findcalls
        hyperparams = struct('batName',batName,'Date', saveDate,'micType', micType, 'callNum',...
        callNum,'recNum',recNum, 'saveTime',saveTime, 'callType',callType, 'lowpass', lowpass,'highpass',highpass,'sessionType',...
        sessionType,'Shift',shift,'sessID',sessID);
        findcalls(senv,wd,fs,findcall_inputs{:});
        % hyperparams should be appended onto struct
        % run each cut through calltrigger.
        append_parameters(wd,hyperparams,FullStamps)
    end


    end
end
 function append_parameters(wd,hyperparams,triggers)
 iterative_dir = strcat(wd,'Analyzed_auto\');
 iterative_length = dir([iterative_dir,'*',Subj,'*',Date,'*',ExpTime,'*.mat']);
        
 for xis = 1:length(iterative_length)
     iterative_name=iterative_length(xis).name;
     retreat = load([iterative_dir iterative_name],'callpos');
     %retreat = matfile([iterative_dir iterative_name],'Writable', true);
     cutcall = retreat.callpos;
     istrigger = trigger_or_not(cutcall,triggers);
     %retreat.istrigger = istrigger;
     %retreat.hyperparams = hyperparams;
     clear retreat;
     save([iterative_dir iterative_name],'istrigger','hyperparams','-append', '-v6')
 end    
 end
function data_raw = load_audio_data(wd,filename,fileType,dataVarName)

switch fileType
    case 'wav'
        data_raw = audioread(fullfile(wd, filename));
    case 'mat'
        data_raw = load(fullfile(wd, filename));
        if isfield(data_raw,'analyzed')
            if data_raw.analyzed
                data_raw = [];
                return
            end
        end
        data_raw = data_raw.(dataVarName);
end
end
end



function [istrigger] = trigger_or_not(callpos,triggers)
if isempty(triggers)
    istrigger=0;
    return
end
length_triggers = length(triggers);
delim = '_';
trigger_stamp= [];
for ix =1:length_triggers
    trigger_stamp(ix,:) = strsplit(string(triggers{ix,1}),delim);
end
upper_range = callpos(2); % set boundaries to parse between
lower_range = callpos(1);
% check and see if callpos is within first or last index of triggers
for ss = 1:length_triggers
    extract_stamp = trigger_stamp(ss,2);
    elementary = intersect(find(lower_range >= extract_stamp),find(extract_stamp<=upper_range));

end
if isempty(elementary)
    istrigger = 0;
else
    istrigger = 1 ;
end
end


%%% Added Functions for Tobias call-cutting
function [date_array,name_array,time_array] = unique_gen(wd)
check_files = dir(fullfile(wd, '*.wav'));
% this generates a list of all the wav files


names = [check_files.name];
string_holder = strsplit(names,'.wav');
string_holder = string_holder';
% now i have Nx1 cell of vocs, now get rid of everything after the
% truedelim.
%length_names = length(names);
% split each string by whitespace

% split into before and after good delim
%extract delimiter based on set position _4
% strpos is 17, so cut everything after 17
%delim_true = '_VocTrigger';
%string_holder_true = [];

for i = 1:length(string_holder)
    if strlength(string_holder(i,1)) < 17
      % do nothing (don't take in the entry since it is empty)  
    else
    string_holder_true(i) = extractBefore(string(string_holder(i,1)),17);
    end
end

unique_array = unique(string_holder_true);
unique_array(1) = [];
unique_array = unique_array';
for is = 1:length(unique_array)
    date_array(is) = extractBetween(unique_array(is,1),6,11);
    name_array(is) = extractBetween(unique_array(is,1),1,4);
    time_array(is) = extractBetween(unique_array(is,1),13,16);
end
% now find unique values inside each array
% return date_array, name_array, and time_array.

% now we need to check and create the unique date,time and subj
% combinations

end
