% MCS raw data processing. in progress... 
% check
% https://github.com/multichannelsystems/McsMatlabDataTools/wiki/Tutorial
% for info on MCS raw data loading into MATLAB
% andrea3vs@gmail.com
% Sep 2023 14
% routine 1. loads the .h5 data, creates a channel map, plots first min of
% raw data, hipass filters, converts to int16 and saves for kilosort

% next steps: 1. create channel map only once! 2. automatize Kilosort 3.
% refine sorting/clustering parameters, 4.
% find a way to extract the waveforms
% (https://github.com/cortex-lab/KiloSort/issues/35)

clc; clear all; close all;

%% check if MCS packages are installed

funName = 'McsData.m'; 

if exist(funName, 'file') ~= 2
    disp('Function McsData.m not found. Downloading from GitHub...');    
    url = 'https://github.com/multichannelsystems/McsMatlabDataTools/archive/refs/heads/master.zip';
    prevDir = cd(userpath);
    zip_file = websave('MCSpackage.zip', url);
    unzip(zip_file, userpath);
    addpath(genpath(fullfile(userpath,'McsMatlabDataTools-master')));
    cd(prevDir)
    delete(zip_file);
    fprintf('MCS package downloaded into %s and added to MATLAB path\n',userpath)
end

%% grab data

% load .h5 data
[datafile,datadir] = uigetfile({'*.h5'},'select .h5 file','C:\Users\User\Documents');
% datafile        = '2023-05-10T15-06-45McsRecording-BTC3.h5';
fileFullPath    = fullfile(datadir,datafile);
filename        = strsplit(datafile,'.h5');
filename        = filename{1};
oldir           = cd(datadir);
rawData         = McsHDF5.McsData(fileFullPath);                            % load raw data, requires MATLAB McsMatlabDataTools package
day             = datestr(now, 'yyyymmdd_HHMM');
subdir4KS       = [datadir day '_4KS3'];
filenOut        = [filename '_raw01' '.bin'];
mkdir(subdir4KS);
% prompt 

% Define the prompt
prompt = {'Load all recording? (''y'' or range in sec, e.g. 0-10):',...
    'Enter batch size in sec (if previous is a range, this won''t matter) :'};
dlgtitle = 'Input';
dims = [1 45];
definput = {'y','5'};
response = inputdlg(prompt, dlgtitle, dims, definput);
batchSz = str2double(response{2}); % batch size in sec
if isempty(response)
    disp('No input provided. Defaulting to loading all recording');
    response2 = 'y';
else
    response2 = response{1};
end
if strcmpi(response2, 'y')
    disp('Loading all recording...');
    allRecTime = true;
else
    % Check if the response is a number
    numSec = str2double(strsplit(response2,'-'));
    if isnan(numSec)
        disp('Invalid input. Please enter ''y'' or a range of sec (two numbers separated by an hyphen).');
    else
        fprintf('Loading only %.1f-%.1f seconds\n', numSec);
        allRecTime = false;
    end
end

perc        = 10;
    
% channel map
name        = '256MEA100-30iR-ITO_kilosortChanMap';
Nchannels   = 256;
connected   = true(Nchannels, 1);
% chanMap     = (1:Nchannels)';
% chanMap0ind = chanMap - 1;
chars           = 'A':'R';
chars           = chars(chars~='I'& chars~='Q');                            % there's no I or Q in the electrodes ID...
chars           = cellstr(chars(:));
numbrs          = 1:16;
numbrs          = arrayfun(@num2str,numbrs,'UniformOutput',false);
[letters, numbers] = meshgrid(chars,numbrs);
layoutSites     = strcat(letters(:), numbers(:)); 
siteDistX       = 100;                                                      % in microns
siteDistY       = 100;
spacingX        = 0:siteDistX:siteDistX*(numel(numbrs)-1);
spacingY        = 0:siteDistY:siteDistY*(numel(numbrs)-1);
[Xcoord, Ycoord] = meshgrid(spacingX,spacingY);
channelMap1     = table(layoutSites(:),num2cell(Xcoord(:)+1),num2cell(Ycoord(:)+1),'VariableNames',...
    {'channelID','X_coord','Y_coord'});
xcoords         = Xcoord(:)+1;
ycoords         = Ycoord(:)+1;
kcoords         = ones(Nchannels,1);                                        % grouping of channels (i.e. tetrode groups)

% check the DataSubType attribute of the streams
for stream = rawData.Recording{1}.AnalogStream
    sub_type = stream{1}.DataSubType;
    if strcmp(sub_type, 'Electrode')
        electrode_data = stream{1};
    elseif strcmp(sub_type, 'Auxiliary')
        analog_data = stream{1};
    elseif strcmp(sub_type, 'Digital')
        digital_data = stream{1};
    end
end
% check data
% figure; plot(rawData.Recording{1}.AnalogStream{1},[])

% channel_index is the row index in the ChannelData matrix
% channel_index   = rawData.Recording{1}.AnalogStream{1}.Info.ChannelID;
% group_id        = rawData.Recording{1}.AnalogStream{1}.Info.GroupID(channel_index);

% % 96 well plates, I assume?
% rows_96 = reshape(repmat(('A':'H')',1,12)',1,96);
% wells_96 = arrayfun(@(x)({[rows_96(x) num2str(mod(x-1,12)+1)]}), 1:96);
% well = wells_96(group_id + 1);

% position        = rawData.Recording{1}.AnalogStream{1}.Info.Label(channel_index); % Channel Position in the Well
position        = rawData.Recording{1}.AnalogStream{1}.Info.Label;              % grab position of the channels in the well plate

if size(position,1) == 1
    load(fullfile('D:\MEA_data','Default256MEA_map.mat'))
end
chanMap         = zeros(Nchannels,1);
nonConnected    ={'A1','R1','A16','R16'};                                   % these 4 channels are used as internal reference
nonConnectedNum = Nchannels-numel(nonConnected)+1:1:Nchannels;              % assign last 4 positions (dataStream will have only 1:252)
n = 1;
for j = 1:Nchannels
    index0      = find(matches(position,layoutSites{j}));
    if isempty(index0)
        chanMap(j)  = nonConnectedNum(n);
        n = n+1;
    else
        chanMap(j)  = index0;                                               % 253:256 need to be added as empty channals to the data before writing to bin file
    end      
end
chanMap0ind     = chanMap -1;
connected       = ~any(chanMap == nonConnectedNum,2);
channelMap2     = table(chanMap,chanMap0ind,connected,'VariableNames',{'chanMap','chanMap0ind','Connected'});
channelMap1     = [channelMap1,channelMap2];
writetable(channelMap1,[name '.xlsx'])
sampleTimeStamps = rawData.Recording{1}.AnalogStream{1}.ChannelDataTimeStamps;
tick            = rawData.Recording{1}.AnalogStream{1}.Info.Tick(1);            % time difference between two samples in microsec
recPhaseStarts  = [1 find(diff(sampleTimeStamps) > tick)];
deltaT          = mean(diff(sampleTimeStamps))./1e6;                            % deltaT in sec
fs              = 1/deltaT;                                                     % acquisition freq

% Save channel map as .mat

save(fullfile(subdir4KS, [name '.mat']), ...
    'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'name', 'fs')

% Save as json for KS4

chanMapStruct = struct('chanMap',chanMap0ind,'xc',xcoords,'yc',ycoords,'kcoords',kcoords,'n_chan',numel(chanMap0ind));
jsonChanMap = jsonencode(chanMapStruct);
jsonName    = fullfile(subdir4KS, [name '.json']);
fileID      = fopen(jsonName, 'w');
fprintf(fileID, '%s', jsonChanMap);
fclose(fileID);

% plot raw data. maybe change to just few min as sanity check?
fileSz      = size(electrode_data.ChannelDataTimeStamps,2);
totSec      = fileSz/fs; 
nSamBatch   = floor(batchSz*fs);
nBatch      = floor(fileSz/nSamBatch);
segLength   = totSec / nBatch;
BatchSamp   = cell(1, nBatch);
% generate batch boundaries in sec
for i = 1:nBatch
    startVal    = (i - 1) * segLength;
    endVal      = i * segLength;
    BatchSamp{i} = [startVal, endVal];
end

if ~allRecTime
    BatchSamp = {numSec}; % override batch size if a range is specified
end

f1g = figure('Name',[filename, '_qualityCheck'],'Color','w','NumberTitle','off','WindowState','maximized');
tiledlayout('flow');

%% load data in batches

InfoOut = struct('chanMap',channelMap1,'acqFq',fs);
progBar = waitbar(0, 'Processing data...');
for j = 1:numel(BatchSamp)

    cfg.window      = BatchSamp{j};  % Specify the time window in seconds
    partialStream   = electrode_data.readPartialChannelData(cfg);
    if j == numel(BatchSamp)
        indxT           = 1:numel(partialStream.ChannelDataTimeStamps); % read the last timepoint only if it's the last batch (otherwise it would be repeated)
    else
        indxT           = 1:numel(partialStream.ChannelDataTimeStamps)-1;
    end

    rawData1            = partialStream.ChannelData(:,indxT);
    timeVector          = double(partialStream.ChannelDataTimeStamps(indxT))./1e6;
    % to millvolts
    if all(contains(electrode_data.DataUnit,'pV'))
        rawData1        = rawData1.*1E-09;
    else
        sprintf('cannot scale to mV, data is %s. not pV!',electrode_data.DataUnit{1,1})
    end
    % timeVector      = 0:deltaT:(size(rawData1,2)-1)*deltaT;

    baselineMin = 2; % baseline in min
    if baselineMin*60 > batchSz
        baselineMin = batchSz/60;
    end
    sampIndx = floor(fs*(baselineMin*60));
    
    % check if the batch is smaller than entire rec?
    if sampIndx > size(rawData1,2)
        sampIndx = size(rawData1,2);
    end

    if j==1
        meanVolt        = mean(rawData1(:,sampIndx),'all');
        ChanlFlactuations = rms(rawData1(:,1:sampIndx)');
        chanRMS = NaN(size(connected));
        chanRMS(connected) = ChanlFlactuations;
        chanRMS2d       = reshape(chanRMS,numel(unique(xcoords)),[]);
        nexttile
        imagesc(chanRMS2d);
        axis('square')
        c1 = colorbar;
        c1.Label.String = 'RMS (mV)';
        maxNoise        = prctile(ChanlFlactuations,[0,99.5],'all');
        corruptChan     = ChanlFlactuations>maxNoise(2);

    end

    % CAR
    avgRaw1     = median(rawData1,1,'omitnan');
    rawData2    = rawData1 - avgRaw1;
    rawData2(corruptChan,:) = 0;            % corrupted channels to zero volts
    % low pass at 8K?
    HPFrq       = 8e3; 
    if fs> HPFrq*2
        rawData2 = Digitalfilter(rawData2, fs, 0.4e4 ,'TYPE', 'low', 'ORDER', 2);
    end

    % prepare data to be saved as int16
    maxFqPS = 1e4;
    pspectr1 = NaN(sum(connected),round(fs/2.5)+1);
    timeRangePS     = 10;                                                       % first 30s as baseline for PS?
    if size(rawData1,2)/fs < timeRangePS
        timeRangePS = floor(size(rawData1,2)/fs);
    end
    numSamplePS     = floor(fs*timeRangePS);
    for nch = 1:sum(connected)        
        yTmp        = rawData1(nch,1:numSamplePS);
        [pspectr1(nch,:),frq]       = pspectrum(yTmp,fs,'power','FrequencyResolution',10,'FrequencyLimits',[0,maxFqPS]);
    end
    PSall = 10*log10(pspectr1);
    PSall(isinf(PSall)) = NaN;

    % plot all channels PS
    nexttile
    imagesc('CData',PSall,'XData',frq);
    c2 = colorbar;
    c2.Label.String = 'dB';
    ylabel('Channel')
    xlabel('Frequency (Hz)')
    xlim([0 max(frq)])
    ylim([1,sum(connected) ])
    title(sprintf('powerspectrum %.1f-%.1f min',BatchSamp{j}/60));
    caxis([-90, -30])

    % plot PS average +-std
    medianPS = median(PSall,1,'omitnan');
    stdevPS     = std(PSall,[],1,'omitnan');
    curve1 = medianPS + stdevPS;
    curve2 = medianPS - stdevPS;
    x2 = [frq', flip(frq')];  % Concatenate x and its reverse
    inBetween = [curve1, fliplr(curve2)];
    nexttile
    fill(x2, inBetween, [.7,.7,.7],'EdgeColor','none');  % Light green shading
    hold on;  % Keep the plot active
    plot(frq, medianPS, 'k', 'LineWidth', 1);  % Plot the mean curve in red
    xlabel('Frq (Hz)');  % Add labels if needed
    ylabel('dB');
    title(sprintf('powerspectrum %.1f-%.1f min',BatchSamp{j}/60));

    % add empty channels. use them to save average MUA instead of zeros
    % Loop through the padIndices and insert new rows at the specified locations
    for i = 1:numel(nonConnectedNum)
        idx = nonConnectedNum(i);

        % Shift the existing rows down to make room for the new row
        rawData2(idx+1:end+1, :) = rawData2(idx:end, :);

        % Insert the new row
        rawData2(idx, :) = avgRaw1;
    end
    clear rawData1
    
    % if all rec > prepare data for KS, if not just plot 
    if allRecTime
        [MUAforKS3,nanidx,Gain] = sing2int16(rawData2,ChanlFlactuations,perc);                    % convert to 1nt16

        fidOut      = fopen(fullfile(subdir4KS,filenOut), 'a'); % append
        fwrite(fidOut,MUAforKS3(:),'*int16');
        fclose(fidOut);
        InfoOut(j).timeVec  = timeVector;
        InfoOut(j).RMS      = ChanlFlactuations;
        InfoOut(j).Gain2mV  = Gain;
    end
    waitbar(j/numel(BatchSamp), progBar, sprintf('Processed batch %d/%d', j, numel(BatchSamp)));
end
if allRecTime
    save(fullfile(subdir4KS,[filename '_Info.mat']),'InfoOut');
end
close(progBar)
cd(oldir)
saveas(f1g,fullfile(subdir4KS,[f1g.Name, '.fig']))
saveas(f1g,fullfile(subdir4KS,[f1g.Name, '.jpg']))
close(f1g)
appFile = fullfile(which('MCS_analysis_may2024.m'));
resulz  = compiler.build.standaloneApplication(appFile);

% compiler.build.standaloneApplication('MCS_analysis_may2024.m', ...
% 'ExecutableName','MCS_h52bin')

%% functions

function data = Digitalfilter(data, fs, FC, varargin)
%TDTDIGITALFILTER  applies a digital filter to continuous data
%   data = TDTdigitalfilter(DATA, STREAM, FC, 'parameter', 'value', ... ),
%   where DATA is the output of TDTbin2mat, STREAM is the name of the
%   stream store to filter, FC is the cutoff frequency. If FC is a two-
%   element vector, a bandpass filter is applied by default.
%
%   data    contains updated STREAM data store with digital filter applied
%
%   'parameter', value pairs
%       'TYPE'      string, specifies the TYPE of filter to use
%                       'band': bandpass filter (default if FC is two element)
%                       'stop': bandstop filter
%                       'low': low pass butterworth (default if FC is scalar)
%                       'high': high pass butterworth 
%       'ORDER'     scalar, filter order for high pass and low pass filters
%                       (default = 2). If the high pass cutoff frequency is
%                       low enough (usually below ~2 Hz), a custom filter
%                       is applied.
%       'NOTCH'     scalar or array, frequencies to apply notch filter. The
%                       notch filters are always 1st order filters.
%
%   Example
%      data = TDTbin2mat('C:\TDT\OpenEx\Tanks\DEMOTANK2\Block-1');
%      data = TDTdigitalfilter(data, 'Wav1', [300 5000], 'NOTCH', [60 120]);
%      data = TDTdigitalfilter(data, 'Wav2', 10, 'TYPE', 'high', 'ORDER', 4);
%      data = TDTdigitalfilter(data, 'Wav3', 'NOTCH', 60);
%

% poolobj = gcp('nocreate'); % get the pool object, and do it avoiding creating a new one.
% if isempty(poolobj) % check if there is not a pool.
%     c = parcluster;
%     poolsize = c.NumWorkers;
%     parpool(poolsize, 'IdleTimeout',Inf); % create a new pool with the previous poolsize and new specs.
% else
%     % poolsize = poolobj.NumWorkers;
%     % delete( gcp('nocreate')); % delete the current pool object.
% end

% defaults
TYPE = 'band';
ORDER = 2;
NOTCH = [];

% catch if no FC is used, only notch
if strcmp(FC, 'NOTCH')
    varargin = [{'NOTCH'}, varargin];
    TYPE = 'NULL';
    FC = [];
end

% parse varargin
VALID_PARS = {'TYPE','ORDER','NOTCH'};
for ii = 1:2:length(varargin)
    if ~ismember(upper(varargin{ii}), VALID_PARS)
        error('%s is not a valid parameter. See help TDTdigitalfilter', upper(varargin{ii}));
    end
    eval([upper(varargin{ii}) '=varargin{ii+1};']);
end

% validate inputs
if length(FC) == 1
    if ~strcmp(TYPE, 'high') && ~strcmp(TYPE, 'low')
        warning('invalid TYPE for scalar FC, assuming ''low''')
        TYPE = 'low';
    end
elseif length(FC) == 2
    if ~strcmp(TYPE, 'band') && ~strcmp(TYPE, 'stop')
        warning('invalid TYPE for two-dimensional vector, assuming ''band''')
        TYPE = 'band';
    end
end

if ~isempty(NOTCH)
    for i = 1:numel(NOTCH)
        BW = 0.05;
        data = Digitalfilter(data, fs, NOTCH(i)*[1-BW,1+BW], 'TYPE', 'stop', 'ORDER', 1);
        
    end
end

if strcmp(TYPE, 'band')
    data = Digitalfilter(data, fs, FC(2), 'TYPE', 'low', 'ORDER', ORDER);
    data = Digitalfilter(data, fs, FC(1), 'TYPE', 'high', 'ORDER', ORDER);
    return
end

if strcmp(TYPE, 'NULL')
    return
end

Alpha = single(1-exp(-6.283 * FC(1) / 24414.0625));
if strcmp(TYPE, 'high') && Alpha <= 0.0005
    if FC(1) > 0
        %%% Emulate MCSmooth HP filter
        % fprintf('Using alpha smoothing for the high pass filter\n');
        for chan = 1:size(data, 1)
            if size(data,1) == 1
                r2 = zeros(size(data));
                r2(1) = data(1);
            else
                r2 = zeros(size(data(chan,:)));
                r2(1) = data(chan,1);
            end
            for i = 2:length(data)
                r2(i) = Alpha*data(chan,i) + (1-Alpha)*r2(i-1);
            end
            data(chan,:) = data(chan,:) - r2;
        end
    end
else
    Fs = fs; %sampling rate
    [Z, P, K] = butter(ORDER, FC./(Fs/2), TYPE);    
    SOS = zp2sos(Z, P, K);
    data = double(data);
    % use filtfilt here to remove phase distortion
    data = sosfilt(SOS,data,2);
end
end

function varargout = sing2int16(dataSingle,rmschan,perc)
% function converts single data to Int16, to keep resolution,
% discards values > perc*95th percentile of the rmschan (rms
% calc for all individual channels), re-scales to the full
% range of int16 the max and min of the data.

I16r    = single([-2^15,2^15-1]);               % int16 range
avRms   = perc*prctile(rmschan,95);                                        % 10 x 95th percentile of rms distribution
% clip data by discarding all val >rmsMUAm, convert to int16                                  % copy
dataSingle(dataSingle>avRms|dataSingle<-avRms) = NaN;                          % set val that exceed 8 x 95-th perc. of rms to zero
dataSingle = dataSingle.*1e3;                                                   % convert to microvolts to reduce the float points
% stretch the range before converting to int16 to keep resolution!
maxData     = max(dataSingle(:));
minData     = min(dataSingle(:));
% % formula1: new_val = (val - min_val)./ (max_val - min_val).* (new_max - new_min) + new_min  >   non-uniform
% dataSingle  = (dataSingle-minData)./(maxData-minData).*(max(I16r)-min(I16r)) + min(I16r);
% formula2: scale_factor = (new_max-new_min) / (max_val-min_val)     >   uniform
gainFactor  = (max(I16r)-min(I16r))./(maxData-minData);                     % gain is used to stretch uniformly the data
dataSingle  = (dataSingle-minData) .* gainFactor + min(I16r);
% remove Nans
varargout{2} = isnan(dataSingle);               % where the data was clipped
varargout{3} = gainFactor;                      % gain, NB data is not in microV x gain!
dataSingle(isnan(dataSingle)) = 0;
% convert to int16 and apply channel map
varargout{1} = int16(dataSingle);

end