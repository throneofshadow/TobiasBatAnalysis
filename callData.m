classdef callData
    % CallData Version 2.0 Maimon Rose 06/30/17
    % Includes the following experiment types: Maimon's ephys, Deafened bats,
    % Lesioned bats, and Tobias' automatic training.
    %
    % Construct a callData object as follows:
    %
    % cData = callData(fs,expType)
    % fs - audio sampling rate in Hz
    % expType - string correspoding to experiment type
    
    properties(SetAccess = public)
        baseDirs % directories to search for call files
        batNums % bat ID numbers (as strings)
        dateFormat % format for importing dates as datetimes
        birthDates % DOBs for all bats
        treatment
        treatmentType
        lesionedBats % which bats are lesioned (logical)
        nBats % number of bats
        cepstralFlag = false % perform cepstral analysis?
        spCorrFlag = false % perform spCorr analysis?
        batName %name of bat pairs in training
        batchNum % for experiments in separate batches
        callNum %call number (subset of recNum)
        recNum %recording number for the training session
        callTrigger %dt for autoTrain
        callType %'callOnly' or 'callTrigger'
        sessionType %experiment type in the autoTrain
        sessionID %numerical counter for each sessionTye
        micType %which box/microphone used
        xrun %dropped sample number
        callID %unique ID for each call
        callEcho % flag to load either calls or echolocations
        loggerNum % audio logger serial number
        manual_call_class % variable to store potenial manual classification
        lowPass
        highPass
        trig
        shift
        fs % sampling rate
        expType % String indicating which experiment we are dealing with
        
        callWF % call waveform
        expDay % experimental date (in datetime format)
        recTime % time of recording (in datetime format)
        daysOld % number of days from birth for this call
        batNum % bat ID number correspoding to this call
        fName % file path to file from which call was cut
        fName_cut % file path to file with cut call
        nCalls % total number of calls
        callLength % length of call
        file_call_pos % % (for 'ephys') position in samples during session
        callPos % (for 'ephys') position in samples during session
        yinF0 % 'best' fundamental freq. calculated using the yin algorithm
        spCorrF0 % mean fundamental freq. calculated using the spCorr algorithm
        weinerEntropy % mean WE
        spectralEntropy % mean SE
        centroid % mean SE
        energyEntropy % mean EE
        RMS % mean RMS
        ap0 % 'best' coarse aperiodicity
        pitchGoodness
        cepstralF0 % fundamental freq. calculated using a cepstral algorithm
        
        % parameters as above, broken down into windows during the call:
        yinF0Win
        ap0Win
        weinerEntropyWin
        spectralEntropyWin
        centroidWin
        energyEntropyWin
        RMSWin
        pitchGoodnessWin
        cepstralF0Win
        spCorrF0Win
        
        maxCalls = 150e3
        minF0 = 200 % in Hz
        maxF0 = 20e3 % in Hz
        
        windowLength = 1.5e-3 % in ms
        integrationWindow = 1.5e-3
        overlap = 1e-3 % in ms
        thresh = 0.01
        yinF0_interp = false
        
        loadWF = true % Load in all call waveforms or not
    end
    properties (Dependent)
        oldestAge
        youngestAge
    end
    methods
        function cData = callData(fs, expType, varargin) % Function to initialize, load, and perform calculations
            cData.fs = fs;
            cData.expType = expType;
            
            % initialize all parameters:
            
            cData = load_call_data(cData,varargin{:});
            
            [cData.yinF0, cData.spCorrF0,...
                cData.weinerEntropy, cData.spectralEntropy, cData.centroid,...
                cData.energyEntropy, cData.RMS, cData.ap0,cData.pitchGoodness,...
                cData.cepstralF0] = deal(zeros(cData.nCalls,1));
            
            [cData.yinF0Win, cData.ap0Win, cData.weinerEntropyWin,...
                cData.spectralEntropyWin, cData.centroidWin,cData.energyEntropyWin,...
                cData.RMSWin, cData.pitchGoodnessWin,cData.cepstralF0Win]...
                = deal(cell(cData.nCalls,1));
            
            tapers = dpss(cData.windowLength*cData.fs,1.5);
            
            if ~cData.loadWF
                cData.callLength = zeros(cData.nCalls,1);
            end
            
            % iterate through all calls
            lastProgress = 0;
            for call_k = 1:cData.nCalls
                
                if cData.loadWF % if we already loaded all the data into this object
                callWF = cData.callWF{call_k};
                else % if we need to load this data on the fly
                    if ~isempty(cData.fName{call_k})
                        callWF = loadCallWF_onTheFly(cData,call_k);
                        cData.callLength(call_k) = length(callWF)/fs;
                    end
                
                end
                
                % Assemble parameter structure for yin algorithm
                nSamples = length(callWF);
                wSize = round(cData.integrationWindow*cData.fs);
                overlap = round(cData.overlap*cData.fs);
                yinParams = struct('thresh',cData.thresh,'sr',cData.fs,'wsize',wSize,...
                    'hop',wSize - overlap,'range',[1 nSamples],...
                    'bufsize',nSamples+2,'APthresh',2.5,...
                    'maxf0',cData.maxF0,'minf0',cData.minF0);
                
                %addpath('C:\Users\phyllo\Documents\MATLAB\yin\')
                %addpath('C:\Users\Julie\workspace\neurobat-callData')
                addpath('C:\Users\bassp\Downloads\Python Work\Rotation_3\Matlab_Code\yin\yin\')
                [f0, ap0] = calculate_yin_F0(callWF,yinParams,cData.yinF0_interp);
                cData.yinF0Win{call_k} = f0;
                cData.ap0Win{call_k} = ap0;
                
                % take the F0 for the window of the call with the "best"
                % (i.e. lowest) aperiodicity) **No longer using this
                % method**
                 [~, idx] = min(ap0);
                 cData.yinF0(call_k) = f0(idx);
                 cData.ap0(call_k) = ap0(idx);
                cData.yinF0(call_k) = nanmedian(f0);
                cData.ap0(call_k) = nanmedian(ap0);
                
                % Assemble parameter structure for windowed feature
                % calculation
                featuresParams = struct('windowLength',cData.windowLength,'fs',cData.fs,...
                    'overlap',cData.overlap,'cepstralFlag',cData.cepstralFlag,...
                    'spCorrFlag',cData.spCorrFlag,'tapers',tapers,...
                    'maxF0',cData.maxF0,'minF0',cData.minF0);
                
                % calculate those features
                
                [cData.weinerEntropyWin{call_k},cData.spectralEntropyWin{call_k},...
                    cData.centroidWin{call_k},cData.energyEntropyWin{call_k},...
                    cData.RMSWin{call_k},cData.pitchGoodnessWin{call_k},...
                    cData.cepstralF0Win{call_k}, cData.spCorrF0Win{call_k}]...
                    = getCallFeatures(callWF,featuresParams);
                
                % store average values of those features
                cData.weinerEntropy(call_k) = mean(cData.weinerEntropyWin{call_k});
                cData.spectralEntropy(call_k) = mean(cData.spectralEntropyWin{call_k});
                cData.centroid(call_k) = mean(cData.centroidWin{call_k});
                cData.energyEntropy(call_k) = mean(cData.energyEntropyWin{call_k});
                cData.RMS(call_k) = mean(cData.RMSWin{call_k});
                
                if cData.cepstralFlag % if we want to use the cepstral algorithm
                    % **No longer using this method**
                    % [pg, idx] = max(cData.pitchGoodnessWin{call_k});
                    % cData.pitchGoodness(call_k) = pg;
                    % cData.cepstralF0(call_k) = cData.cepstralF0Win{call_k}(idx);
                    cData.pitchGoodness(call_k) = nanmedian(cData.pitchGoodnessWin{call_k});
                    cData.cepstralF0(call_k) = nanmedian(cData.cepstralF0Win{call_k});
                end
                
                if cData.spCorrFlag
                    cData.spCorrF0(call_k) = nanmedian(cData.spCorrF0(call_k));
                end
                
                progress = 100*(call_k/cData.nCalls);
                
                if mod(progress,10) < mod(lastProgress,10)
                    fprintf('%d %% of calls processed\n',round(progress));
                end
                
                lastProgress = progress;
                
            end
            
        end
        function n = numArgumentsFromSubscript(obj,~,~)
            n = numel(obj);
        end
        function cData = load_call_data(cData,varargin)
            
            switch cData.expType
                case 'ephys' % Maimon's juvenile ephys experiments or Wujie's adult ephys experiments
                    
                    if ~isempty(varargin)
                        juv_adult = varargin{1};
                    else
                        juv_adult = 'juvenile';
                    end
                    
                    eData = ephysData(juv_adult);
                    
                    cData.baseDirs = eData.baseDirs;
                    cData.batNums = eData.batNums;
                    cData.dateFormat = eData.dateFormat;
                    cData.birthDates = eData.birthDates;
                    cData.nBats = length(cData.batNums);
                    
                    [cData.callWF, cData.batNum, cData.fName] = deal(cell(cData.maxCalls,1));
                    [cData.daysOld, cData.callLength] = deal(zeros(cData.maxCalls,1));
                    cData.callPos = zeros(cData.maxCalls,2);
                    cData.file_call_pos = zeros(cData.maxCalls,2);
                    
                    cData.expDay = datetime([],[],[]);
                    nlg_rec_str = 'neurologger_recording';
                    
                    if isempty(cData.callEcho)
                        cData.callEcho = input('Call or Echo?','s');
                    end
                    if strcmp(cData.callEcho,'call')
                        cutFName = 'cut_call_data.mat';
                    elseif strcmp(cData.callEcho,'echo')
                        cutFName = 'cut_echo_data.mat';
                    end
                    
                    call_k = 1;
                    for b = 1:cData.nBats % iterate across all experimental bats
                        switch juv_adult
                            case 'juvenile'
                                call_info_str = 'juv_call_info';
                                nlg_dir_str = [cData.baseDirs{b} 'bat' cData.batNums{b} filesep '*neurologger*'];
                            case 'adult'
                                call_info_str = 'call_info';
                                nlg_dir_str = [cData.baseDirs{b} '*neurologger*'];
                        end
                        
                        nlgDirs = dir(nlg_dir_str);
                        for d = 1:length(nlgDirs) % iterate across all recording days
                            audioDir = [fullfile(nlgDirs(d).folder,nlgDirs(d).name) filesep 'audio\ch1\'];
                            
                            if strcmp(cData.callEcho,'echo')
                                call_info_fName = [audioDir call_info_str '_' cData.batNums{b} '_echo.mat'];
                                if ~exist(call_info_fName,'file')
                                    continue
                                end
                                
                                s = load(call_info_fName);
                                call_info = s.call_info;
                                
                                s = load([audioDir cutFName]);
                                cut_call_data = s.cut_call_data;
                                cut_call_data = cut_call_data(~[cut_call_data.noise]);
                                
                                assert(all([cut_call_data.uniqueID] == [call_info.callID]));
                                
                                echo_calls = arrayfun(@(x) strcmp(x.echoCall,'juvEcho'),call_info);
                                cut_call_data = cut_call_data(echo_calls);
                            else
                                % get call data and time in recording
                                % (corrected for clock drift)
                                s = load([audioDir cutFName]);
                                cut_call_data = s.cut_call_data;
                                if isempty(cut_call_data)
                                   continue 
                                end
                                if strcmp(juv_adult,'adult')
                                    
                                    exp_day_idx = strfind(nlgDirs(d).name,nlg_rec_str) + length(nlg_rec_str);
                                    expDate = nlgDirs(d).name(exp_day_idx:exp_day_idx+length(cData.dateFormat)-1);
                                    
                                    call_info_fname = dir([audioDir 'call_info_' cData.batNums{b} '_' cData.callEcho '_' expDate '.mat']);
                                    
                                    if isempty(call_info_fname)
                                       continue 
                                    end
                                    
                                    if length(call_info_fname) > 1
                                        keyboard
                                    end
                                    s = load(fullfile(call_info_fname.folder,call_info_fname.name));
                                    call_info = s.call_info;
                                    
                                    cut_call_data = cut_call_data(~[cut_call_data.noise]);
                                    
                                    batIdx = unique(cellfun(@(call) find(cellfun(@(bNum) strcmp(bNum,cData.batNums{b}),call)),{cut_call_data.batNum}));
                                    
                                    if length(batIdx) == 1
                                        callpos = horzcat(cut_call_data.corrected_callpos);
                                        callpos = callpos(batIdx,:);
                                        [cut_call_data.corrected_callpos] = deal(callpos{:});
                                    else
                                        keyboard
                                    end
                                    
                                    assert(all([cut_call_data.uniqueID] == [call_info.callID]));
                                    
                                    bat_calls = cellfun(@(x) ischar(x{1}) && contains(x,cData.batNums{b}),{call_info.behaviors});
                                    cut_call_data = cut_call_data(bat_calls);
                                    
                                end
                            end
                            
                            if ~isempty([cut_call_data.f_num])
                                
                                cutCalls = {cut_call_data.cut};
                                call_pos_expDay = vertcat(cut_call_data.corrected_callpos)/1e3; % convert to seconds
                                file_call_pos_expDay = vertcat(cut_call_data.callpos);
                                call_length_expDay = cellfun(@length, cutCalls)/cData.fs;
                                exp_day_idx = strfind(nlgDirs(d).name,nlg_rec_str) + length(nlg_rec_str);
                                expDatetime = datetime(nlgDirs(d).name(exp_day_idx:exp_day_idx+length(cData.dateFormat)-1),'InputFormat',cData.dateFormat);
                                
                                for c = 1:length(cutCalls)
                                    if ~cut_call_data(c).noise && ~any(isnan(cut_call_data(c).corrected_callpos))
                                        cData.callWF{call_k} = cutCalls{c};
                                        cData.callLength(call_k) = call_length_expDay(c);
                                        cData.callPos(call_k,:) = call_pos_expDay(c,:);
                                        cData.file_call_pos(call_k,:) = file_call_pos_expDay(c,:);
                                        cData.expDay(call_k,1) = expDatetime;
                                        cData.batNum{call_k} = cData.batNums{b};
                                        cData.daysOld(call_k) = days(expDatetime - cData.birthDates{b});
                                        cData.fName{call_k} = cut_call_data(c).fName;
                                        cData.callID(call_k) = cut_call_data(c).uniqueID;
                                        call_k = call_k + 1;
                                    end
                                end
                            end
                        end
                    end
                    cData.nCalls = call_k-1;
                    
                    
                case 'lesion'
                    
                    maxBatchNums = 5;
                    batchNums = 1:maxBatchNums;
                    cData.loadWF = false;
                    s = load('E:\lesion_recordings\all_lesion_bat_info.mat','all_lesion_bat_info');
                    batInfo = s.all_lesion_bat_info;
                    
                    cData.dateFormat = 'yyyyMMdd';
                    dateRegExp = '_\d{8}T';
                    
                    [cData.callWF, cData.fName, cData.treatment] = deal(cell(cData.maxCalls,1));
                    [cData.callLength, cData.daysOld, cData.batchNum] = deal(zeros(cData.maxCalls,1));
                    cData.expDay = datetime([],[],[]);
                    cData.callID = [1:cData.maxCalls]';
                    
                    treatmentGroups = fieldnames(batInfo)';
                    
                    if ~isempty(varargin)
                        varargin = varargin{1};
                        for subs_k = 1:2:length(varargin)
                            switch varargin{subs_k}
                                
                                case 'treatment'
                                    treatmentGroups = treatmentGroups(ismember(treatmentGroups,varargin{subs_k+1}));
                                case 'batchNum'
                                    batchNums = varargin{subs_k+1};
                            end
                            
                        end
                    end
                    
                    
                    call_k = 1;
                    for t = 1:length(treatmentGroups)
                        groupStr = treatmentGroups{t};
                        batchNumsGroup = find(~structfun(@isempty,batInfo.(groupStr).birthDates))';
                        batchNumsGroup = batchNumsGroup(ismember(batchNumsGroup,batchNums));
                        
                        for b = batchNumsGroup
                            bStr = ['batch' num2str(b)];
                            avg_birth_date = mean(batInfo.(groupStr).birthDates.(bStr));
                            callFiles = dir([batInfo.(groupStr).baseDir  batInfo.(groupStr).treatmentType.(bStr) '*' bStr '*_Call_*.mat']);
                            if ~isempty(callFiles)
                                for c = 1:length(callFiles)
                                    cData.treatment{call_k} = batInfo.(groupStr).treatmentType.(bStr);
                                    exp_date_str = regexp(callFiles(c).name,dateRegExp,'match');
                                    exp_date_str = exp_date_str{1}(2:end-1);
                                    cData.expDay(call_k) = datetime(exp_date_str,'inputFormat',cData.dateFormat);
                                    if isdatetime(avg_birth_date)
                                        cData.daysOld(call_k) = days(cData.expDay(call_k) - avg_birth_date);
                                    end
                                    cData.fName{call_k} = [batInfo.(groupStr).baseDir callFiles(c).name];
                                    cData.batchNum(call_k) = b;
                                    call_k = call_k + 1;
                                end
                            end
                        end
                    end
                    cData.nCalls = call_k-1;
                    
                    
                case 'deafened'
                    maxBatchNums = 5;
                    batchNums = 1:maxBatchNums;
                    cData.loadWF = false;
                    s = load('E:\deafened_recordings\all_deaf_bat_info.mat','all_deaf_bat_info');
                    batInfo = s.all_deaf_bat_info;
                    
                    cData.dateFormat = 'yyyyMMddHHmmss';
                    dateRegExp = '_\d{8}T';
                    timeRegExp = 'T\d{6}_';
                    
                    [cData.callWF, cData.fName, cData.fName_cut, cData.treatment] = deal(cell(cData.maxCalls,1));
                    [cData.callLength, cData.daysOld, cData.batchNum] = deal(zeros(cData.maxCalls,1));
                    cData.expDay = datetime([],[],[]);
                    cData.callID = [1:cData.maxCalls]';
                    cData.file_call_pos = zeros(cData.maxCalls,2);
                    
                    treatmentGroups = fieldnames(batInfo)';
                    
                    if ~isempty(varargin)
                        varargin = varargin{1};
                        for subs_k = 1:2:length(varargin)
                            switch varargin{subs_k}
                                
                                case 'treatment'
                                    treatmentGroups = treatmentGroups(ismember(treatmentGroups,varargin{subs_k+1}));
                                case 'batchNum'
                                    batchNums = varargin{subs_k+1};
                            end
                            
                        end
                    end
                    
                    
                    call_k = 1;
                    for t = 1:length(treatmentGroups)
                        groupStr = treatmentGroups{t};
                        batchNumsGroup = find(~structfun(@isempty,batInfo.(groupStr).birthDates));
                        batchNumsGroup = batchNumsGroup(ismember(batchNumsGroup,batchNums))';
                        group_data_baseDir = batInfo.(groupStr).baseDir(1:strfind(batInfo.deaf.baseDir,'all_cut_calls')-1);
                        
                        for b = batchNumsGroup
                            bStr = ['batch' num2str(b)];
                            avg_birth_date = mean(batInfo.(groupStr).birthDates.(bStr));
                            callFiles = dir([batInfo.(groupStr).baseDir  batInfo.(groupStr).treatmentType.(bStr) '*' bStr '*_Call_*.mat']);
                            raw_rec_dir = [group_data_baseDir 'all_raw_recordings\' batInfo.(groupStr).data_dir_str.(bStr) filesep];
                            if ~isempty(callFiles)
                                for c = 1:length(callFiles)
                                    cData.treatment{call_k} = batInfo.(groupStr).treatmentType.(bStr);
                                    exp_date_str = regexp(callFiles(c).name,dateRegExp,'match');
                                    time_date_str = regexp(callFiles(c).name,timeRegExp,'match');
                                    exp_datetime_str = [exp_date_str{1}(2:end-1) time_date_str{1}(2:end-1)];
                                    cData.expDay(call_k) = datetime(exp_datetime_str,'inputFormat',cData.dateFormat);
                                    if isdatetime(avg_birth_date)
                                        cData.daysOld(call_k) = days(cData.expDay(call_k) - avg_birth_date);
                                    end
                                    cData.fName_cut{call_k} = [batInfo.(groupStr).baseDir callFiles(c).name];
                                    cData.fName{call_k} = [raw_rec_dir callFiles(c).name(1:strfind(callFiles(c).name,'_Call_')-1) '.mat'];
                                    if ~exist(cData.fName{call_k},'file')
                                        cData.fName{call_k} = [];
                                    end
                                    cData.batchNum(call_k) = b;
                                    s = load(cData.fName{call_k},'callpos');
                                    cData.file_call_pos(call_k,:) = s.callpos;
                                    call_k = call_k + 1;
                                end
                            end
                        end
                    end
                    cData.nCalls = call_k-1;
                    cData.expDay = cData.expDay';
                case 'autoTrain'
                    cData.loadWF = false;
                    %cData.baseDirs = 'Y:\users\tobias\vocOperant\box7\bataudio\preShift\Analyzed_auto\';
                    %cData.baseDirs = 'Y:\users\tobias\vocOperant\box7\bataudio\shift_20-25\triggered\Analyzed_auto\';
                    %cData.baseDirs = 'Y:\users\tobias\vocOperant\box1\bataudio\Iterative\';
                    cData.baseDirs = 'C:\Users\bassp\Downloads\Python Work\Rotation 1\reexamples\h5files\';
                    %cData.baseDirs = 'Y:\users\tobias\vocOperant\box7\bataudio\noShift\Detected_calls\Translated_files\Analyzed_auto\';
                    %cData.baseDirs = 'C:\Users\tobias\Desktop\analysis\bataudio\call groups\all\';
                    %cData.baseDirs = 'C:\Users\tobias\Desktop\analysis\bataudio\April2017\cut_preprocessed\';
                    %cData.baseDirs = 'C:\Users\tobias\Desktop\analysis\bataudio\test\';
                    %files = dir(cData.basedir);
                    % Get a logical vector that tells which is a directory.
                    %dirFlags = [files.isdir];
                    % Extract only those that are directories.
                    %subFolders = files(dirFlags);
                    %for k = 1 : length(subFolders) % remember to add an
                    %end, also! add subfolders+Analyzed_auto\ to callFiles flag
                    % add check for Analyzed_auto in dirs as well so we
                    % dont iterate over folders without Analyzed_auto
                    callFiles = dir([cData.baseDirs '*.mat']);
                    % WRITE CASE FOR HYPERPARAMS VS NONHYPERPARAMS,
                    % findFIELD or SEEKSTRUCT
                    
                    cData.nCalls = length(callFiles);
                    [cData.fName, cData.batName, cData.callType, cData.micType, cData.sessionType] = deal(cell(cData.nCalls,1)); % initialize call data cells
                    [cData.recNum, cData.callNum, cData.xrun, cData.sessionID,] = deal(zeros(cData.nCalls,1)); % initialize call data arrays
                    %cData.expDay = datetime([],[],[]);
                    for call_k = 1:length(callFiles)
                        s = load([cData.baseDirs callFiles(call_k).name]);
                        si_x = isfield(s,'hyperparams');
                        if si_x == 1
                            cData.fName{call_k} = callFiles(call_k).name;
                            
                            try
                                cData.batName{call_k} = s.hyperparams.batName;
                            catch
                                cData.batName{call_k} = 'CoEd';
                            end
                            try
                                cData.sessionType{call_k} = s.sessionType;
                            catch
                                cData.sessionType{call_k} = 0;
                            end
                            %cData.callWF{call_k} = s.rawData';
                            cData.callNum(call_k) = s.hyperparams.callNum;
                            
                            
                            cData.lowPass(call_k) = (s.hyperparams.lowpass);
                            try
                                cData.highPass(call_k)= str2double(s.hyperparams.highpass);
                            catch
                                cData.highPass(call_k) = 0;
                            end
                            cData.trig(call_k) = s.istrigger;
                            cData.shift(call_k) = s.hyperparams.Shift;
                            sessID = s.hyperparams.sessID;
                            if ischar(sessID)
                                sessID = str2double(sessID);
                            end
                            if size(sessID) == [1 2]
                                sessID = 0;
                            end
                            cData.sessionID(call_k) = sessID;
                            if isfield(s,'xrun')
                                cData.xrun(call_k) = s.hyperparams.xrun;
                            else
                                cData.xrun(call_k) = 0;
                            end
                            if isfield(s,'callTrigger')
                                cData.exp
                                expDay(call_k) = s.hyperparams.Date;
                            else
                                formatout = 'yyyymmdd';
                                A = datestr(datenum(s.hyperparams.Date,'yymmdd',2000),formatout);
                                %s.Date = datestr(s.Date,dateformat);
                                cData.expDay(call_k) = str2double(A);
                            end
                            try
                                cData.callType{call_k} = s.hyperparams.callType;
                            catch
                                if cData.trig == 1
                                    cData.callType{call_k} = 'VocTrigger';
                                else
                                    cData.callType{call_k} = 'NoTrigger';
                                end
                            end
                            cData.micType{call_k} = s.hyperparams.micType;
                        
                            cData.expDay = cData.expDay';
                        else
                            cData.fName{call_k} = callFiles(call_k).name;
                            
                            try
                                cData.batName{call_k} = s.batName;
                            catch
                                cData.batName{call_k} = 'CoEd';
                            end
                            try
                                cData.sessionType{call_k} = s.sessionType;
                            catch
                                cData.sessionType{call_k} = 0;
                            end
                            %cData.callWF{call_k} = s.rawData';
                            try
                                cData.callNum(call_k) = s.callNum;
                            catch
                                cData.callNum(call_k) = 0;
                            end
                            
                            cData.lowPass(call_k) = (s.lowpass);
                            try
                                cData.highPass(call_k)= str2double(s.highpass);
                            catch
                                cData.highPass(call_k) = 0;
                            end
                            cData.trig(call_k) = s.istrigger;
                            cData.shift(call_k) = s.Shift;
                            sessID = s.sessID;
                            if ischar(sessID)
                                sessID = str2double(sessID);
                            end
                            if size(sessID) == [1 2]
                                sessID = 0;
                            end
                            cData.sessionID(call_k) = sessID;
                            if isfield(s,'xrun')
                                cData.xrun(call_k) = s.xrun;
                            else
                                cData.xrun(call_k) = 0;
                            end
                            if isfield(s,'callTrigger')
                                cData.expDay(call_k) = s.Date;
                            else
                                formatout = 'yyyymmdd';
                                A = datestr(datenum(s.Date,'yymmdd',2000),formatout);
                                %s.Date = datestr(s.Date,dateformat);
                                cData.expDay(call_k) = str2double(A);
                            end
                            try
                                cData.callType{call_k} = s.callType;
                            catch
                                if cData.trig == 1
                                    cData.callType{call_k} = 'VocTrigger';
                                else
                                    cData.callType{call_k} = 'NoTrigger';
                                end
                            end
                            try
                                cData.micType{call_k} = s.micType;
                            catch
                                cData.micType{call_k}=0;
                            end
                        end
                        cData.expDay = cData.expDay';
                    end
                case 'pratData'
                    cData.loadWF = false;
                    cData.baseDirs = 'E:\Yossi_vocalization_data\';
                    
                    [batInfo,batMetadata] = get_yossi_bat_info(cData.baseDirs);
                    
                    treatmentTypes = {'isolated','control'};
                    
                    cData.nCalls = sum(cellfun(@(x) length(vertcat(batInfo.(x)(:).callPos)),treatmentTypes));
                    [cData.fName, cData.treatment, cData.batNum] = deal(cell(cData.nCalls,1));
                    [cData.callLength, cData.daysOld, cData.treatment] = deal(zeros(cData.nCalls,1));
                    cData.recTime = datetime(zeros(0,3));
                    cData.batNums = cellfun(@(x) unique({batInfo.(x).batNum}),treatmentTypes,'un',0);
                    cData.batNums = [cData.batNums{:}];
                    cData.file_call_pos = zeros(cData.nCalls,2);
                    
                    cData.birthDates = datetime(zeros(0,3));
                    for b = 1:length(cData.batNums)
                        dob = batMetadata.DOB{batMetadata.ID==str2double(cData.batNums{b})};
                        cData.birthDates(b) = dob;
                    end
                    cData.callID = [1:cData.nCalls]';
                    call_k = 1;
                    for t = treatmentTypes
                        treat = t{1};
                        for file_k = 1:length(batInfo.(treat))
                            for file_call_k = 1:size(batInfo.(treat)(file_k).callPos,1)
                                cData.treatment{call_k} = treat;
                                cData.fName{call_k} = [cData.baseDirs batInfo.(treat)(file_k).fName];
                                cData.recTime(call_k) = batInfo.(treat)(file_k).recTime;
                                cData.batNum{call_k} = batInfo.(treat)(file_k).batNum;
                                cData.daysOld(call_k) = round(days(batInfo.(treat)(file_k).recTime - cData.birthDates(strcmp(cData.batNums,cData.batNum{call_k}))));
                                cData.treatmentType(call_k) = batInfo.(treat)(file_k).treatmentType;
                                cData.file_call_pos(call_k,:) = batInfo.(treat)(file_k).callPos(file_call_k,:);
                                call_k = call_k + 1;
                            end
                        end
                    end
                    
                case 'piezo_recording'
                    
                    cData.loadWF = false;
                    dataDir = 'C:\Users\phyllo\Documents\Maimon\acoustic_recording\';
                    audio_base_dir = 'Z:\users\Maimon\acoustic_recording\audio\';
                    recordingLogs = readtable([dataDir 'recording_logs.csv']);
                    recordingLogs = recordingLogs(logical(recordingLogs.usable),:);
                    treatmentGroups = readtable([dataDir 'bat_info.csv']);
                    bat_ID_table_idx = contains(recordingLogs.Properties.VariableNames,'Bat_');
                    nExp = size(recordingLogs,1);

                    cData.batNums = cellfun(@num2str,num2cell(treatmentGroups.BatNum),'un',0);
                    
                    cData.nBats = length(cData.batNums);
                    cData.birthDates = cell(1,cData.nBats);
                    
                    for b = 1:cData.nBats
                        dob = treatmentGroups.DOB(b);
                        cData.birthDates{b} = dob;
                    end
                    
                    [cData.callWF, cData.batNum, cData.treatment, cData.fName, cData.fName_cut] = deal(cell(cData.maxCalls,1));
                    [cData.loggerNum, cData.callLength, cData.daysOld] = deal(zeros(cData.maxCalls,1));
                    cData.callPos = zeros(cData.maxCalls,2);
                    cData.file_call_pos = zeros(cData.maxCalls,2);
                    cData.expDay = datetime([],[],[]);
                    
                    date_str_format = 'mmddyyyy';
                    cutFName = 'cut_call_data.mat';
                    analysis_dir_name = 'Analyzed_auto';
                    
                    call_k = 1;
                    
                    for d = 1:nExp % iterate across all recording days
                        expDate = recordingLogs.Date(d);
                        exp_day_bat_nums = recordingLogs{d,bat_ID_table_idx};
                        exp_day_logger_nums = recordingLogs{d,find(bat_ID_table_idx)+1};
                        logger_bat_ID_idx = ~isnan(exp_day_bat_nums) & ~isnan(exp_day_logger_nums) & ~ismember(exp_day_logger_nums,str2double(recordingLogs.malfunction_loggers{1}));
                        
                        exp_day_logger_nums = exp_day_logger_nums(logger_bat_ID_idx);
                        exp_day_bat_nums = exp_day_bat_nums(logger_bat_ID_idx);
                        
                        dateStr = datestr(expDate,date_str_format);
                        audioDir = fullfile(audio_base_dir,dateStr,'audio\ch1\');
                        
                        cut_call_data = load(fullfile(audioDir,cutFName));
                        cut_call_data = cut_call_data.cut_call_data;
                        all_cut_call_files = dir(fullfile(audioDir,analysis_dir_name,'*Call*.mat'));
                        
                        assert(length(all_cut_call_files) == length(cut_call_data));
                        
                        AL_info = load(fullfile(audioDir,'AL_class_info.mat'));
                        
                        cut_call_data = cut_call_data(AL_info.usableIdx);
                        all_cut_call_files = all_cut_call_files(AL_info.usableIdx);

                        predicted_loggerID = predict_AL_identity(AL_info);
                        noise_idx = cellfun(@length,predicted_loggerID) == 0;
                        
                        cut_call_data = cut_call_data(~noise_idx);
                        predicted_loggerID = predicted_loggerID(~noise_idx);
                        all_cut_call_files = all_cut_call_files(~noise_idx);
                        
                        identifiable_call_idx = cellfun(@length,predicted_loggerID) == 1;
                        usable_loggerID = [predicted_loggerID{identifiable_call_idx}];
                        loggerID = nan(1,sum(~noise_idx));
                        loggerID(identifiable_call_idx) = AL_info.logger_nums(usable_loggerID);
                        
                        if isempty(cut_call_data)
                            continue
                        end
                        
                        call_pos_expDay = vertcat(cut_call_data.corrected_callpos)/1e3; % convert to seconds
                        file_call_pos_expDay = vertcat(cut_call_data.callpos);
                        call_length_expDay = arrayfun(@(x) length(x.cut)/cData.fs,cut_call_data);
                        
                        for c = 1:length(cut_call_data)
                            cData.callLength(call_k) = call_length_expDay(c);
                            cData.callPos(call_k,:) = call_pos_expDay(c,:);
                            cData.file_call_pos(call_k,:) = file_call_pos_expDay(c,:);
                            cData.expDay(call_k,1) = expDate;
                            cData.batNum{call_k} = num2str(exp_day_bat_nums(exp_day_logger_nums == loggerID(c)));
                            cData.loggerNum(call_k) = loggerID(c);
                            cData.fName{call_k} = cut_call_data(c).fName;
                            cData.fName_cut{call_k} = fullfile(all_cut_call_files(c).folder,all_cut_call_files(c).name);
                            cData.callID(call_k) = cut_call_data(c).uniqueID;
                            
                            if ~isempty(cData.batNum{call_k})
                                b = strcmp(cData.batNums,cData.batNum{call_k});
                                cData.daysOld(call_k) = days(expDate - cData.birthDates{b});
                                cData.treatment{call_k} = treatmentGroups.Treatment{b};
                            else
                                cData.batNum{call_k} = 'multiple_bats_detected';
                                cData.treatment{call_k} = 'multiple_bats_detected';
                            end
                            call_k = call_k + 1;
                        end
                    end
                    cData.nCalls = call_k-1;
                    
                    
                otherwise
                    ME = MException('CallData:inputError','Experiment Type not recognized');
                    throw(ME);
            end
            
            % if we initialized different data structures to have a length
            % of 'maxCalls' go ahead and shorten those to remove empty
            % elements
            callProperties = properties(cData)';
            for prop = callProperties
                if all(size(cData.(prop{:})) == [cData.maxCalls,1])
                    cData.(prop{:}) = cData.(prop{:})(1:cData.nCalls);
                elseif size(cData.(prop{:}),1) == cData.maxCalls && size(cData.(prop{:}),1) > 1
                    cData.(prop{:}) = cData.(prop{:})(1:cData.nCalls,:);
                end
            end
            
%            cData.callID = reshape(cData.callID,[1 cData.nCalls]);
            
        end
        function varargout = subsref(cData,S)
            if length(S) == 2
                switch S(1).type
                    case '()'
                        nSubs = length(S(1).subs);
                        if ~rem(nSubs,2)
                            callIdx = true(cData.nCalls,1);
                            for idx = 1:2:nSubs
                                switch S(1).subs{idx}
                                    case 'cellInfo'
                                        callIdx = callIdx & strcmp(cData.cellInfo,S(1).subs{idx+1});
                                    case 'daysOld'
                                        daysOldIdx = false(cData.nCalls,1); %make an index assuming all 0 initially
                                        for d = S(1).subs{idx+1} %for the input that you're searching for
                                            daysOldIdx = daysOldIdx | round(cData.daysOld)==d; %make true if the day of indexed call is listed
                                        end
                                        callIdx = callIdx & daysOldIdx;
                                    case 'expDay'
                                        if length(S(1).subs{idx+1}) == 1
                                            callIdx = callIdx & (cData.expDay == S(1).subs{idx+1});
                                        else
                                            callIdx = callIdx & (cData.expDay >= S(1).subs{idx+1}(1) & cData.expDay < S(1).subs{idx+1}(2));
                                        end
                                    case 'batNum'
                                        batNumIdx = false(cData.nCalls,1);
                                        for b = S(1).subs{idx+1}
                                            batNumIdx = batNumIdx | strcmp(cData.batNum,b);
                                        end
                                        callIdx = callIdx & batNumIdx;
                                    case 'treatment'
                                        treatmentIdx = false(cData.nCalls,1);
                                        for b = S(1).subs{idx+1}
                                            treatmentIdx = treatmentIdx | strcmp(cData.treatment,b);
                                        end
                                        callIdx = callIdx & treatmentIdx;
                                    case 'batName'
                                        batNameIdx = false(cData.nCalls,1);
                                        for bN = S(1).subs{idx+1}
                                            batNameIdx = batNameIdx | ~cellfun(@isempty, strfind(cData.batName,bN)); %cellfun(@(x) ~isempty(strfind(x,bN)),cData.batName);
                                        end
                                        callIdx = callIdx & batNameIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.batName);
                                    case 'batchNum'
                                        batchNumIdx = false(cData.nCalls,1);
                                        for bN = S(1).subs{idx+1}
                                            batchNumIdx = batchNumIdx | cData.batchNum == bN;
                                        end
                                        callIdx = callIdx & batchNumIdx;
                                    case 'callID'
                                        callIdx = callIdx & ismember(cData.callID,S(1).subs{idx+1});
                                    case 'callNum'
                                        callNumIdx = false(cData.nCalls,1);
                                        for rN = S(1).subs{idx+1}
                                            callNumIdx = callNumIdx | cData.callNum==rN;
                                        end
                                        callIdx = callIdx & callNumIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.callNum);
                                    case 'callTrigger'
                                        callIdx = callIdx & (cData.callTrigger >= S(1).subs{idx+1}(1) & cData.callTrigger <= S(1).subs{idx+1}(2));
                                    case 'callType'
                                        callTypeIdx = false(cData.nCalls,1);
                                        for cT = S(1).subs{idx+1}
                                            callTypeIdx = callTypeIdx | cellfun(@(x) ~isempty(strfind(x,cT)),cData.callType);
                                        end
                                        callIdx = callIdx & callTypeIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.callType);
                                    case 'micType'
                                        callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.micType);
                                    case 'recNum'
                                        recNumIdx = false(cData.nCalls,1);
                                        for rN = S(1).subs{idx+1}
                                            recNumIdx = recNumIdx | cData.recNum==rN;
                                        end
                                        callIdx = callIdx & recNumIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.recNum);
                                    case 'sessionID'
                                        sessionIDIdx = false(cData.nCalls,1);
                                        for sI = S(1).subs{idx+1}
                                            sessionIDIdx = sessionIDIdx | cData.sessionID==sI;
                                        end
                                        callIdx = callIdx & sessionIDIdx;
                                        %                                        callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.sessionID);
                                    case 'callLength'
                                        callLengthIdx = false(cData.nCalls,1);
                                        for cL = S(1).subs{idx+1}
                                            callLengthIdx = callLengthIdx | cData.callLength==cL;
                                        end
                                        callIdx = callIdx & callLengthIdx;
                                    case 'sessionType'
                                        sessionTypeIdx = false(cData.nCalls,1);
                                        for sT = S(1).subs{idx+1}
                                            sessionTypeIdx = sessionTypeIdx | strcmp(sT,cData.sessionType);
                                        end
                                        callIdx = callIdx & sessionTypeIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.sessionType);
                                    case 'lowPass'
                                        lowPassIdx = false(cData.nCalls,1);
                                        for sT = S(1).subs{idx+1}
                                            lowPassIdx = lowPassIdx | cData.lowPass==sT;
                                        end
                                        callIdx = callIdx & lowPassIdx;
                                    case 'highPass'
                                        highPassIdx = false(cData.nCalls,1);
                                        for sT = S(1).subs{idx+1}
                                            highPassIdx = highPassIdx | cData.highPass == sT;
                                        end
                                        callIdx = callIdx & highPassIdx;
                                    case 'trig'
                                        trigIdx = false(cData.nCalls,1);
                                        for sT = S(1).subs(idx+1)
                                            trigIdx = trigIdx | cData.trig == sT;
                                        end
                                        callIdx = callIdx & trigIdx;
                                    case 'shift'
                                        shiftIdx = false(cData.nCalls,1);
                                        for sT = S(1).subs(idx+1)
                                            shiftIdx = shiftIdx | cData.shift == sT;
                                        end
                                        callIdx = callIdx & shiftIdx
                                    case 'xrun'
                                        callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.xrun);
                                    otherwise
                                        display('indexing variable not recognized')
                                        return
                                end
                            end
                            if iscell(cData.(S(2).subs)(callIdx)) && ~any(cellfun(@ischar,cData.(S(2).subs)(callIdx)))
                                try
                                    varargout = {vertcat(cData.(S(2).subs){callIdx})};
                                catch
                                    try
                                        varargout = {[cData.(S(2).subs){callIdx}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                                
                            else
                                varargout = {cData.(S(2).subs)(callIdx,:)};
                            end
                        else
                            display('Indexing in VocalData must come in pairs');
                            return
                        end
                    otherwise
                        switch S(2).type
                            case '{}'
                                try
                                    varargout = {vertcat(cData.(S(1).subs){S(2).subs{:}})};
                                catch
                                    try
                                        varargout = {[cData.(S(1).subs){S(2).subs{:}}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                            otherwise
                                try
                                    varargout = {builtin('subsref',cData,S)};
                                catch err
                                    switch err.message
                                        case 'Too many output arguments.'
                                            builtin('subsref',cData,S);
                                        otherwise
                                            display(err)
                                            return
                                    end
                                end
                        end
                end
            else
                try
                    varargout = {builtin('subsref',cData,S)};
                catch err
                    switch err.message
                        case 'Too many output arguments.'
                            builtin('subsref',cData,S);
                        otherwise
                            display(err)
                            return
                    end
                end
            end
            
            
        end
        function [interCallInterval, interDayIdx] = getICI(cData)
            cPos = cData.callPos;
            interCallInterval = [0; cPos(2:end,1) - cPos(1:end-1,2)]; % ICI(k) = time between call(k) and previous call
            interDayIdx = abs(diff([cData.expDay(1); cData.expDay]))>duration(1,0,0) | interCallInterval<0;
        end
        function [boutDuration, boutICI] = get_bout_duration(cData)
            boutSeparation = 1;
            
            cLength = cData.callLength;
            
            [interCallInterval, interDayIdx] = getICI(cData);
            interBoutIdx = [find(interCallInterval>boutSeparation); length(interCallInterval)];
            
            boutDuration = zeros(1,length(interBoutIdx)-1);
            boutICI = cell(1,length(interBoutIdx)-1);
            
            for bout_k = 1:length(interBoutIdx)-1
                last_call_in_bout = interBoutIdx(bout_k+1)-1;
                calls_in_bout = interBoutIdx(bout_k)+1:last_call_in_bout;
                if length(calls_in_bout) > 2
                    if any(interDayIdx(calls_in_bout))
                        last_call_in_bout = calls_in_bout(find(interDayIdx(calls_in_bout),1,'first'))-1;
                    end
                    calls_in_bout = interBoutIdx(bout_k)+1:last_call_in_bout;
                    boutDuration(bout_k) = sum(cLength(calls_in_bout)) + sum(interCallInterval(calls_in_bout));
                    boutICI{bout_k} = interCallInterval(calls_in_bout);
                    if boutDuration(bout_k) < 0
                        keyboard
                    end
                else
                    boutDuration(bout_k) = NaN;
                end
                
            end
            boutICI = vertcat(boutICI{:});
            boutDuration = boutDuration(boutDuration~=0 & ~isnan(boutDuration));
        end
        function oldestAge = get.oldestAge(obj)
            oldestAge = max(obj.daysOld);
        end
        function youngestAge = get.youngestAge(obj)
            youngestAge = min(obj.daysOld);
        end
        function callWF = getCallWF(cData,call_k)
            if cData.loadWF
                callWF = cData.callWF{call_k};
            else
                callWF = loadCallWF_onTheFly(cData,call_k);
            end
        end
        function plotSpectrogram(cData,axisHandle,specParams,callInfo)
            
            if length(callInfo) == 1
                if cData.loadWF
                    data = cData.callWF{callInfo};
                else
                    data = loadCallWF_onTheFly(cData,callInfo);
                end
            else
                data = callInfo;
            end
            
            specWin = kaiser(cData.fs*cData.windowLength,0.5);
            nfft = 256;
            freqRange = [0 40e3];
            specFreqs = linspace(freqRange(1),freqRange(2),nfft);
            spectrogram(data,specWin,cData.fs*cData.overlap,specFreqs,cData.fs,'yaxis');
            
            addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
            if ~isfield(specParams,'spec_caxis_factor')
                specParams.spec_caxis_factor = 0.75;
            end
            
            if length(data) <= specParams.spec_win_size
                return
            end
            
            [~,f,t,ps] = spectrogram(data,gausswin(specParams.spec_win_size),specParams.spec_overlap_size,specParams.spec_nfft,specParams.fs,'yaxis');
            if length(t) <= 1
                return
            end
            f = f*1e-3; % frequency in kHz
            t = t*1e3;
            t = [0 t(end-1)];
            ps = 10*log10(ps);
            ps(isinf(ps)) = NaN;
            imagesc(axisHandle,t,f,ps)
            cmap = spec_cmap();
            colormap(cmap);
            ylim(specParams.spec_ylims);
            caxis([min(ps(:))*specParams.spec_caxis_factor max(ps(:))]);
            set(axisHandle,'YDir','normal')
            
        end
        function cData = manual_classify_call_type(cData,varargin)
            randomizeFlag = false;
            overWriteFlag = true;
            orig_rec_plot_win = 5;
            specParams = struct('spec_win_size',512,'spec_overlap_size',500,'spec_nfft',1024,'fs',cData.fs,'spec_ylims',[0 60],'spec_caxis_factor',0.75);
            if isempty(cData.manual_call_class)
                cData.manual_call_class = cell(1,cData.nCalls);
            end
            
            if ~isempty(varargin)
                start_call_k = varargin{1};
            else
                start_call_k = 1;
            end
            
            last_orig_wav_data = '';
            call_idx = start_call_k:cData.nCalls;
            
            if ~overWriteFlag
                call_idx = setdiff(call_idx,find(~cellfun(@isempty,cData.manual_call_class)));
            end
            
            if randomizeFlag
                call_idx = call_idx(randperm(length(call_idx)));
            end
            
            for call_k = call_idx
                data = getCallWF(cData,call_k);
                playback_fs = min(cData.fs,200e3);
                sound(data,playback_fs);
                
                origRec_fName = cData.fName{call_k};
                load_orig_wav_data = ~strcmp(origRec_fName,last_orig_wav_data);
                
                if load_orig_wav_data
                    dataFull = audioread(origRec_fName);
                    last_orig_wav_data = origRec_fName;
                end
                
                subplot(2,1,1)
                cla
                plotSpectrogram(cData,gca,specParams,data)
                subplot(2,1,2);
                hold on
                if load_orig_wav_data
                    cla
                    plot((1:length(dataFull))/cData.fs,dataFull,'k');
                end
                plot((cData.file_call_pos(call_k,1)+(0:length(data)-1))/cData.fs,data);
                
                if length(dataFull)/cData.fs > orig_rec_plot_win
                    xlim([cData.file_call_pos(call_k,1)/cData.fs - orig_rec_plot_win/2 cData.file_call_pos(call_k,1)/cData.fs + orig_rec_plot_win/2])
                end
                
                display([datestr(cData.expDay(call_k)) ' Call ' num2str(call_k)]);
                
                repeat = 1;
                repeat_k = 1;
                while repeat
                    class = input('Call type?','s');
                    if isempty(class)
                        pause(0.1);
                        if repeat_k < 3
                            sound(data,playback_fs/(2*repeat_k));
                        else
                            startIdx = max(1,cData.file_call_pos(call_k,1) - (orig_rec_plot_win/2)*cData.fs);
                            endIdx = min(length(dataFull),cData.file_call_pos(call_k,1) + (orig_rec_plot_win/2)*cData.fs);
                            sound(dataFull(startIdx:endIdx),playback_fs);
                            repeat_k = 1;
                        end
                        repeat_k = repeat_k + 1;
                        
                    else
                        if strcmp(class,'n')
                            cData.manual_call_class{call_k} = 'noise';
                            repeat = 0;
                        elseif strcmp(class,'m')
                            cData.manual_call_class{call_k} = 'mating';
                            repeat = 0;
                        elseif strcmp(class,'t')
                            cData.manual_call_class{call_k} = 'trill';
                            repeat = 0;
                        elseif strcmp(class,'o')
                            cData.manual_call_class{call_k} = 'other';
                            repeat = 0;
                        elseif strcmp(class,'i')
                            cData.manual_call_class{call_k} = 'interesting';
                            repeat = 0;
                        elseif strcmp(class,'stop')
                            return
                        elseif strcmp(class,'pause')
                            keyboard
                        end
                    end
                end
            end
            
        end
    end
end

function callWF = loadCallWF_onTheFly(cData,call_k)

switch cData.expType
    case 'pratData'
        callWF = audioread(cData.fName{call_k});
        callWF = callWF(cData.file_call_pos(call_k,1):cData.file_call_pos(call_k,2));
    case 'autoTrain'
        if ~isempty(cData.fName{call_k})
            callWF = load([cData.baseDirs cData.fName{call_k}]);
            % 3 switch statements for 'recbuf','cut', and 'dataRaw'
            % use an ifexist clause
            
            callWF = callWF.cut;
        end
    case 'deafened'
        callWF = load(cData.fName{call_k});
        try
            callWF = callWF.cut;
        catch
            callWF = callWF.finalcut;
        end
        callWF = reshape(callWF,numel(callWF),1);
        
    case 'lesion'
        callWF = load(cData.fName{call_k});
        callWF = callWF.cut;
        callWF = reshape(callWF,numel(callWF),1);
    case 'piezo_recording'
        callWF = load(cData.fName_cut{call_k});
        callWF = callWF.cut;
        callWF = reshape(callWF,numel(callWF),1);        
        
    otherwise
        display('No functionality to load callWF for the experiment type');
        keyboard;
end

end

function [F0,ap] = calculate_yin_F0(callWF,P,interpFlag)

% Adapted from Vidush Mukund vidush_mukund@berkeley.edu
% March 13, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builds acts as a wrapper to call the YIN Algorithm package;
% based on the script written of Yosef Prat
%   P.minf0:    Hz - minimum expected F0 (default: 30 Hz)
%   P.maxf0:    Hz - maximum expected F0 (default: SR/(4*dsratio))
%   P.thresh:   threshold (default: 0.1)
%   P.relfag:   if ~0, thresh is relative to min of difference function (default: 1)
%   P.hop:      samples - interval between estimates (default: 32)
%   P.range:    samples - range of samples ([start stop]) to process
%   P.bufsize:  samples - size of computation buffer (default: 10000)
%	P.sr:		Hz - sampling rate (usually taken from file header)
%	P.wsize:	samples - integration window size (defaut: SR/minf0)
%	P.lpf:		Hz - intial low-pass filtering (default: SR/4)
%	P.shift		0: shift symmetric, 1: shift right, -1: shift left (default: 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% call the yin function from the yin package
nSamples = length(callWF);
R = yin(callWF,P);

% F0 = 2.^(R.f0(round(P.wsize/2/P.hop)+1:end) + log2(440)); % convert from octaves to Hz
% ap = R.ap0(round(P.wsize/2/P.hop)+1:end);
F0 = 2.^(R.f0 + log2(440));
ap = R.ap0;

if interpFlag
    
    ap(ap ~= real(ap)) = abs(ap(ap ~= real(ap)));
    ap = smooth(ap,5);
    
    F0(ap > P.APthresh) = NaN;
    
    % any frequency values that are below the accepted range of fundamentals or
    % above the maximum fundamental frequency are replced with NaN
    F0(F0<=P.minf0 | F0>=P.maxf0 | F0<= 110) = NaN;
    
    % drop all NaN values (i.e. f0 values outside of the accepted range of
    % possible frequency values)
    ap = ap(~isnan(F0));
    F0 = F0(~isnan(F0));
    
    if nargin < 5
        N = length(F0);
    end
    
    if length(F0) > 1
        % perform a 1-D interpolation of the fundamenetal frequencies
        t=(1:length(F0)).*P.hop/P.sr + 0.5/P.minf0;
        tt=linspace(t(1),nSamples/P.sr,2*N+1);
        tt = tt(2:2:end-1);
        F0 = interp1(t,F0,tt);
        ap = interp1(t,ap,tt);
    elseif isempty(F0)
        F0 = NaN;
        ap = NaN;
    end
    idx = ~isnan(F0);
    F0 = F0(idx);
    ap = ap(idx);
end

if isempty(F0)
    F0 = NaN;
    ap = NaN;
end

F0 = smooth(F0,2);
ap = smooth(ap,2);

end
function F0 = calculate_spCorr_F0(callWF,p)
maxLag = round((1/p.minF0)*p.fs);
r = xcorr(callWF, maxLag, 'coeff');
F0 = spPitchCorr(r, p.fs, p.maxF0, p.minF0);
end
function [f0] = spPitchCorr(r, fs, mxf, mnf)
% NAME
%   spPitchCorr - Pitch Estimation via Auto-correlation Method
% SYNOPSIS
%   [f0] = spPitchCorr(r, fs)
% DESCRIPTION
%   Estimate pitch frequencies via Cepstral method
% INPUTS
%   r        (vector) of size (maxlag*2+1)x1 which contains Corr coefficients.
%             Use spCorr.m
%   fs       (scalar) the sampling frequency of the original signal
% OUTPUTS
%   f0       (scalar) the estimated pitch
% AUTHOR
%   Naotoshi Seo, April 2008
% SEE ALSO
%   spCorr.m
% search for maximum  between 2ms (=500Hz) and 20ms (=50Hz)
%adjusted parameters are 4000 Hz and 500 Hz
ms2=floor(fs/mxf); % 2ms
ms20=floor(fs/mnf); % 20ms
% half is just mirror for real signal
r = r(floor(length(r)/2):end);
[maxi,idx]=max(r(ms2:ms20));
f0 = fs/(ms2+idx-1);
end
function [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
    pitchGoodness, cepstralF0, spCorrF0Win] = getCallFeatures(callWF,P)

L = length(callWF);
L_frame = P.windowLength*P.fs;
L_step = L_frame - P.overlap*P.fs;
nFrame = floor((L-(L_frame-L_step))/L_step);

if nFrame > 0
    [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
        pitchGoodness, cepstralF0, spCorrF0Win] = deal(zeros(1,nFrame));
    for fr = 1:nFrame
        frameIdx = ((fr-1)*L_step +1):(((fr-1)*L_step)+L_frame);
        frame = callWF(frameIdx);
        [weinerEntropy(fr), spectralEntropy(fr), centroid(fr),...
            energyEntropy(fr), RMS(fr), pitchGoodness(fr), cepstralF0(fr),...
            spCorrF0Win(fr)] = getFeatures(frame,P);
    end
else
    [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
        pitchGoodness, cepstralF0, spCorrF0Win] = deal(NaN);
end


end
function [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
    pitchGoodness, cepstralF0, spCorrF0] = getFeatures(frame,P)

fs = P.fs;
tapers = P.tapers;
minF0 = P.minF0;
maxF0 = P.maxF0;

[AFFT, F] = afft(frame,fs);

%wiener entropy
weinerEntropy = exp(sum(log(AFFT)) ./ length(AFFT)) ./ mean(AFFT);

%spectral entropy
spectralEntropy = sentropy(AFFT);

%center of spectral mass (centroid)
centroid = sum( F' .* AFFT ) ./ sum(AFFT);

%Energy entropy
energyEntropy = sentropy(frame);

%RMS
RMS = sqrt(mean(frame.^2));

% Cepstral F0 pitch goodness
if P.cepstralFlag
    [pitchGoodness, cepstralF0] = cepstrumAnalysis(frame,tapers,fs);
else
    pitchGoodness = NaN;
    cepstralF0 = NaN;
end

if P.spCorrFlag % if we want to use the spCorr algorithm
    spCorrParams = struct('fs',fs,'minF0',minF0,'maxF0',maxF0);
    spCorrF0 = calculate_spCorr_F0(frame,spCorrParams);
else
    spCorrF0 = NaN;
end


end
function [AFFT, F] = afft(sig,fs,nfft)

if size(sig,1) == 1
    sig = sig';
end

L = size(sig,1);

if nargin < 3 || isempty(nfft)
    nfft = 2^nextpow2(L);
end

F = fs/2*linspace(0,1,nfft/2+1);
AFFT = fft(sig,nfft)./L;
AFFT = 2*abs(AFFT(1:nfft/2+1,:));

end
function ent = sentropy(v)

if size(v,2) == 1
    v = v';
end

v = v + abs(min(min(v(:)),0));
v = v./(ones(size(v,1),1)*sum(v));
v(v==0) = 1;
ent = -sum(v.*log(v));

end
function [pitchGoodness, f0] = cepstrumAnalysis(frame,tapers,fs)

frame = frame.*tapers(:,1);
tmp=rceps(frame);
[pitchGoodness, idx] = max(tmp(25:floor(end/2)));
f0 = (1/idx)*fs;

end

