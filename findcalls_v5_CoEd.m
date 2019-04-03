%

function findcalls_v5_CoEd(wd,fs,varargin,dataVarName)
%% INPUTS: wd refers to working directory
%%         fs refers to audio sampling rate
%%         varargin, dataVarName refers to .mat or not etc
%%         


%ask if wav or mat file and change filetype
if nargin == 2
    wav_mat_file = input('wav (1) or mat (2) file?');
    if wav_mat_file == 1
        fileType = 'wav';
    elseif wav_mat_file == 2
        fileType = 'mat';
    else
        display('Invalid input');
        return
    end
elseif nargin == 3
    fileType = varargin{1};
end

%dataVarName = 'recbuf';
%debug?
debug = 0;
merge = 1;

%call cut parameters
calllength=0.01*fs;%parameter has to be adapted to sig/noise ratio premerge call duration min
%for call eval
durthresh=30; %in ms : postmerge is the call long enough
mergethresh = 20e-3; % in s : any windows less than mergethresh (5ms) apart will be combined
rmsthresh=0.0005;%has to be adapted to recording quality : is it loud enough premerge?
powerRatioThresh = Inf;
high_filter_cutoff = 2000; %high pass filter to eliminate low level noise
low_cenv_filter_cutoff = 500; %smoothies out the envelope


adaptiveThreshold = false;
%changes threshold adaptively based on std of a reshaped data
if adaptiveThreshold
    n_std_thresh = 15;
else
    thresh = 0.005; %determines the initial premerge voltage cutoff min
end

rec_files = dir([wd '*.' fileType]);
%make new directory of analyzed files
anal_dir = [wd 'Analyzed_auto' filesep];
if ~isdir(anal_dir)
    mkdir(anal_dir);
end
n_files = length(rec_files); %loop through each file
%to help with debugging
if debug
    figure; %start image
    manual_classification = zeros(1,n_files); %make variable of each call
    low2high = nan(1,n_files); %range of power for each call
end

total_callcount=0;
for fln = 1:n_files
    file_callcount = 0;
    filename=rec_files(fln).name;
    disp(filename)
    disp(['Analyzing file: ' num2str(fln) ' of ' num2str(n_files)])
    disp('...')
    
    %pull in audio data
    switch fileType
        case 'wav'
            data_raw = audioread([wd filename]);
        case 'mat'
            inputStruct = load([wd filename]);
            data_raw = inputStruct.(dataVarName);
            data_raw = data_raw(:,1);
    end
    
    %high pass to eliminate low noise
    [high_b, high_a] = butter(2,2*high_filter_cutoff/fs,'high'); 
    data = filtfilt(high_b,high_a,data_raw);
    %reshape the data based on std of rounded data
    if adaptiveThreshold
        data_round10 = 10*floor(length(data)/10);
        reshape_data_MA = reshape(data(1:data_round10),[],10);
        thresh = n_std_thresh*min(std(reshape_data_MA));
    end
    
    hilbenv = abs(hilbert(data));
    [b,a] = butter(2,2*low_cenv_filter_cutoff/fs);
    senv = filtfilt(b,a,hilbenv);
    if length(unique(sign(senv-thresh))) > 1
        thresh_up = find(diff(sign(senv-thresh))==2);
        thresh_down = find(diff(sign(senv-thresh))==-2);
        if isempty(thresh_up) || isempty(thresh_down) % check that threhold crossing happens within file and not at beginning or end
            wins = [];
            indx = [];
        else
            if length(thresh_up)==length(thresh_down) && thresh_up(1)<thresh_down(1);
                wins = [thresh_up,thresh_down];
            elseif length(thresh_up)>length(thresh_down)
                wins = [thresh_up(1:end-1),thresh_down];
            elseif length(thresh_down)>length(thresh_up)
                wins = [thresh_up,thresh_down(2:end)];
            elseif length(thresh_up)==length(thresh_down) && thresh_up(1)>thresh_down(1)
                wins = [thresh_up(1:end-1),thresh_down(2:end)];
            end
            if ~isempty(wins)
                win_cell = mat2cell(wins,ones(1,size(wins,1)),2);
                diffs = cellfun(@diff,win_cell);
                rms_win = cellfun(@(x) rms(data(x(1):x(2))),win_cell);
                indx=find((diffs>calllength) & (rms_win>rmsthresh));
            else
                indx = [];
            end
        end
    else
        wins = [];
        indx = [];
    end
    
    if debug
        cla
        hold on
        plot(data);
        plot(wins(indx,:)',max(data)*ones(2,size(indx,1)));
    end
    if ~isempty(indx)
        if merge
            wins = merge_wins(wins(indx,:),fs,mergethresh);
        end
        for w = 1:size(wins,1)
            callpos = wins(w,:);
            cut = data_raw(callpos(1):callpos(2));
            t=(length(cut)/fs)*1000;
            H=rms(cut);
            %             peak2med = max(hilbenv(win(1):win(2)))/median(hilbenv(win(1):win(2)));
            powerRatio = bandpower(cut,fs,[0 5e3])/bandpower(cut,fs,[5e3 10e3]);
            if (t>=durthresh && H>rmsthresh && powerRatio<powerRatioThresh)
                
                total_callcount=total_callcount+1;
                disp(['Call count: ' num2str(total_callcount)])
                if debug
                    plot(callpos',max(data)*ones(2,1),'LineWidth',5);
                    sound(cut,min(fs,200e3));
                    plot(senv);
                    plot(get(gca,'xlim'),[thresh thresh],'k')
                    low2high(fln) = bandpower(cut,fs,[0 5e3])/bandpower(cut,fs,[5e3 10e3]);
                    display(sprintf('ratio of power in 0-5kHz to 10-15kHz = %3.2s',low2high(fln)));
                    keyboard;
                else
                    inputStruct.cut = cut;
                    inputStruct.callpos = callpos;
                    save([anal_dir filename(1:end-4) '_Call_' sprintf('%03d',file_callcount) '.mat'],'-struct','inputStruct');
                end
                file_callcount=file_callcount+1;
                
            end
        end
    end
end

end