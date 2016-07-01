function vis_stream(varargin)
% Display an LSL stream.
%
% Keyboard shortcuts:
%   [up arrow]   : increase the y scale of the time series
%   [down arrow] : decrease the y scale of the time series
%   [right arrow]: increase the displayed time range
%   [left arrow] : decrease the displayed time range
%   [page up]    : go up by one page of channels
%   [page down]  : go down by one page of channels
%
% In:
%   StreamName : Stream to display. The name of the stream that you would like to display.
%
%   TimeScale : Initial time scale in seconds. The time range of the display window;
%               can be changed with keyboard shortcuts (see help). Default=5
%
%   DataScale : Initial scale of the data. The scale of the data, in units between horizontal linesCamera;
%               can be changed with keyboard shortcuts (see help). Default=150
%
%   ChannelRange : Channels to display. The channel range to display. Default=[1:32]
%
%   SamplingRate : Sampling rate for display. This is the sampling rate that is used for plotting, in Hz;
%                  for faster drawing. Default=100
%
%   RefreshRate : Refresh rate for dvisisplay. This is the rate at which the graphics are updated, in Hz.
%                 Default=10
%
%   Rereference : Apply common-average re-referencing to the data. Useful for noisy EEG streams.
%                 Default=false
%
%   PageOffset : Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.
%                Default=0
%
%   Position : Figure position. Allows to script the position at which the figures should appear.
%              This is a 4-element vector of the form [X-offset,Y-offset,Width,Height]
%              with all values in pixels.
%              Default=[]
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-07-10
%
%                                uses portions of vis_dataStreamViewer
%                                (c) 2012 by Tim Mullen

% make sure that everything is on the path and LSL is loaded
if ~exist('arg_define','file')
    addpath(genpath(fileparts(mfilename('fullpath')))); end
if ~exist('env_translatepath','file')
    % standalone case
    lib = lsl_loadlib();
else
    % if we're within BCILAB we want to make sure that the library is also found if the toolbox is compiled
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
end

% handle input arguments
streamnames = unique(find_streams);

% init shared handles and clear existing data in base workspace
evalin('caller',['clear ','lsl_*']) 
[chunkname,...
    chunktime, ...
    buffername,...
    inlet,...
    axCamera,...
    linesCamera,...
    axMOBIlab,...
    linesMOBIlab,...
    opts{2}] = deal({});

% Get entire duration of recording
for k = 1 : length(streamnames)
    if strcmp(streamnames{k},'Length')
        streamnames(k) = [] ;
        break;
    end
    if strcmp(streamnames{k},'Markers')
        streamnames(k) = [] ;
        break;
    end
end

% init Camera and MOBIlab properties
for i = 1 : size(streamnames,2)
    switch streamnames{i}
        case 'Camera'
            opts{i} = arg_define(varargin, ...
                arg({'streamname','StreamName'},streamnames{i},streamnames,'LSL stream that should be displayed. The name of the stream that you would like to display.'), ...
                arg({'bufferrange','BufferRange'},10,[],'Maximum time range to buffer. Imposes an upper limit on what can be displayed.'), ...
                arg({'timerange','TimeRange'},5,[],'Initial time range in seconds. The time range of the display window; can be changed with keyboard shortcuts (see help).'), ...
                arg({'datascale','DataScale'},150,[],'Initial scale of the data. The scale of the data, in units between horizontal linesCamera; can be changed with keyboard shortcuts (see help).'), ...
                arg({'channelrange','ChannelRange'},1:6,[],'Channels to display. The channel range to display.'), ...
                arg({'samplingrate','SamplingRate'},64,[],'Sampling rate for display. This is the sampling rate that is used for plotting; for faster drawing.'), ...
                arg({'refreshrate','RefreshRate'},10,[],'Refresh rate for display. This is the rate at which the graphics are updated.'), ...
                arg({'reref','Rereference'},false,[],'Common average reference. Enable this to view the data with a common average reference filter applied.'), ...
                arg({'standardize','Standardize'},false,[],'Standardize data.'), ...
                arg_nogui({'subplot','Subplot'},[1 2 3 4 5 6],[],'Subplot #, if any.'), ...
                arg_nogui({'pageoffset','PageOffset'},0,[],'Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.'), ...
                arg_nogui({'position','Position'},[],[],'Figure position. Allows to script the position at which the figures should appear.'));
            
        case 'MOBIlab'
            opts{i} = arg_define(varargin, ...
                arg({'streamname','StreamName'},streamnames{i},streamnames,'LSL stream that should be displayed. The name of the stream that you would like to display.'), ...
                arg({'bufferrange','BufferRange'},10,[],'Maximum time range to buffer. Imposes an upper limit on what can be displayed.'), ...
                arg({'timerange','TimeRange'},5,[],'Initial time range in seconds. The time range of the display window; can be changed with keyboard shortcuts (see help).'), ...
                arg({'datascale','DataScale'},150,[],'Initial scale of the data. The scale of the data, in units between horizontal linesCamera; can be changed with keyboard shortcuts (see help).'), ...
                arg({'channelrange','ChannelRange'},1:2,[],'Channels to display. The channel range to display.'), ...
                arg({'samplingrate','SamplingRate'},256,[],'Sampling rate for display. This is the sampling rate that is used for plotting; for faster drawing.'), ...
                arg({'refreshrate','RefreshRate'},10,[],'Refresh rate for display. This is the rate at which the graphics are updated.'), ...
                arg({'reref','Rereference'},false,[],'Common average reference. Enable this to view the data with a common average reference filter applied.'), ...
                arg({'standardize','Standardize'},false,[],'Standardize data.'), ...
                arg_nogui({'subplot','Subplot'},[1 2],[],'Subplot #, if any.'), ...
                arg_nogui({'pageoffset','PageOffset'},0,[],'Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.'), ...
                arg_nogui({'position','Position'},[],[],'Figure position. Allows to script the position at which the figures should appear.'));
    end
end

for DeviceNumber = 1 : size(streamnames,2)
    
    % fix up some arguments
    opts{DeviceNumber}.bufferrange = max(opts{DeviceNumber}.bufferrange,opts{DeviceNumber}.timerange);
    
    % variable names to hold the stream's data (in the base workspace)
    chunkname{DeviceNumber}  = genvarname(['lsl_' opts{DeviceNumber}.streamname '_chunk']);
    buffername{DeviceNumber} = genvarname(['lsl_' opts{DeviceNumber}.streamname '_stream']);
    chunktime{DeviceNumber}  = genvarname(['lsl_' opts{DeviceNumber}.streamname '_timestamps']);
    
    % create a stream inlet
    inlet{DeviceNumber} = create_inlet(opts{DeviceNumber});

    % create the stream data structure in the base workspace
    assignin('base', buffername{DeviceNumber}, create_streambuffer(opts{DeviceNumber}, inlet{DeviceNumber}.info()));
    
    switch streamnames{DeviceNumber}
        case 'Camera'
            % create the figure
            create_figure_Camera(opts{DeviceNumber});
            
            % set up a timer that reads from LSL
            timer_obj_Camera = timer('Period', 1.0/opts{1}.refreshrate,'ExecutionMode','fixedRate','TimerFcn',@on_timer_Camera,'StartDelay',0.2,'Tag',['lsl_' genvarname(opts{1}.streamname) '_timer']);
            start(timer_obj_Camera);
            
        case 'MOBIlab'
            if isempty(opts{2})
                opts{1,2} = opts{1,1};
                chunkname{1,2} = chunkname{1,1};
                chunktime{1,2} = chunktime{1,1};
                buffername{1,2} = buffername{1,1};
                inlet{1,2} = inlet{1,1};
            end
            
            % create the figure
            create_figure_MOBIlab(opts{DeviceNumber});
            
            % set up a timer that reads from LSL
            timer_obj_MOBIlab = timer('Period', 1.0/opts{2}.refreshrate,'ExecutionMode','fixedRate','TimerFcn',@on_timer_MOBIlab,'StartDelay',0.2,'Tag',['lsl_' genvarname(opts{2}.streamname) '_timer']);
            start(timer_obj_MOBIlab);
    end
    
end

% === nested functions (sharing some handles with each other) ===

    % create a new figure and axes
    function create_figure_Camera(opts)
        figure('Tag',['Fig' buffername{DeviceNumber}],'Name',['LSL:Stream''' opts.streamname ''''],'CloseRequestFcn','delete(gcbf)', ...
            'KeyPressFcn',@(varargin)on_key_Camera(varargin{2}.Key));
        for k = 1 : length(opts.subplot)
            axCamera{k} = subplot(length(opts.subplot),1,opts.subplot(k));
            set(gca,'tag',['Camera ' num2str(k)]);
        end
    end

    function create_figure_MOBIlab(opts)
        figure('Tag',['Fig' buffername{DeviceNumber}],'Name',['LSL:Stream''' opts.streamname ''''],'CloseRequestFcn','delete(gcbf)', ...
            'KeyPressFcn',@(varargin)on_key_MOBIlab(varargin{2}.Key));
        for k = 1 : length(opts.subplot)
            axMOBIlab{k} = subplot(length(opts.subplot),1,opts.subplot(k));
        end
    end

    function on_timer_Camera(varargin)
        try 
            % check if the buffer is still there
            if evalin('base',['exist(''' buffername{1} ''',''var'')'])
                
                % === update buffer contents (happens in the base workspace) ===
                
                % pull a new chunk from LSL
                [chunkdata,timestamps] = inlet{1}.pull_chunk();
                assignin('base',chunkname{1},chunkdata);
                assignin('base',chunktime{1},timestamps);            
               
                % append it to the stream buffer
                evalin('base',['[' buffername{1} '.smax,' buffername{1} '.data(:,1+mod(' buffername{1} '.smax:' buffername{1} '.smax+size(' chunkname{1} ',2)-1,' buffername{1} '.pnts))] = deal(' buffername{1} '.smax + size(' chunkname{1} ',2),' chunkname{1} ');']);
                
                % get the updated stream buffer
                stream = evalin('base',buffername{1});
                

                % reformat the stream buffer to contain only the current block that should be displayed
                samples_to_get = min(stream.pnts, round(stream.srate*stream.opts.timerange));
                stream.data = stream.data(:, 1+mod(stream.smax-samples_to_get:stream.smax-1,stream.pnts));
                [stream.nbchan,stream.pnts,stream.trials] = size(stream.data);
                if ~isempty(timestamps)
                    stream.xmax = max(timestamps)-lsl_local_clock(lib);
                elseif ~isfield(stream,'xmax')
                    stream.xmax = 0;
                end
                stream.xmin = stream.xmax - (stream.pnts-1)/stream.srate;
                
                % === data post-processing for plotting ===
                
                % determine channels and samples to display
                plotchans = stream.opts.channelrange + stream.opts.pageoffset*length(stream.opts.channelrange);
                if isempty(plotchans)
                    plotchans = 1:stream.nbchan;
                else
                    plotchans = intersect(1:stream.nbchan,plotchans);
                end
                plotdata = stream.data(plotchans, round(1 : stream.srate/stream.opts.samplingrate : end));
                plottime = linspace(stream.xmin,stream.xmax,size(plotdata,2));
                
                % re-reference
                if stream.opts.reref
                    plotdata = bsxfun(@minus,plotdata,mean(plotdata)); end
                if stream.opts.standardize
                    plotdata = bsxfun(@times,plotdata,1./std(plotdata')'); end
                
                % zero-mean
				tmpdata = plotdata; tmpdata(isnan(tmpdata(:))) = 0;
                plotdata = bsxfun(@minus, plotdata, mean(tmpdata,2));
                
                % arrange for plotting
                plotoffsets = (0:size(plotdata,1)-1)*stream.opts.datascale;
                plotdata = bsxfun(@plus, plotdata', (plotoffsets(end)/2));
                
                
                % === actual drawing ===
                
                % draw the block contents...
                if ~isempty(plotdata)
                    if ~exist('linesCamera','var') || isempty(linesCamera)                        
                        for k = 1 : length(opts{1}.subplot)
                            linesCamera{k} = plot(axCamera{k},plottime,plotdata(:,k));                        
                            ylabel(axCamera{k},'Degrees','FontSize',12);
                            grid(axCamera{k},'minor');
                        end
                        title(axCamera{1},opts{1}.streamname);
                        xlabel(axCamera{length(opts{1}.subplot)},'Time (sec)','FontSize',12);
                    else
                        for k = 1 : length(linesCamera)
                            set(linesCamera{k},'Ydata',plotdata(:,k));
                            set(linesCamera{k},'Xdata',plottime);
                        end
                    end
                
                    % update the axis limit and tickmarks
                    for k = 1 : length(linesCamera)
                        axis(axCamera{k},[stream.xmin stream.xmax -stream.opts.datascale size(plotdata,2)*stream.opts.datascale + stream.opts.datascale]);
                        set(axCamera{k}, 'YTick',(plotoffsets(end)/2), 'YTickLabel', stream.chanlocs(k).labels);
                    end
                end
                
                drawnow;
                
            else
                try 
                    disp(['Deleting timer ' get(timer_obj_Camera,'Tag') '.']);
                catch 
                    disp('Deleting timer.');
                end
                % delete the timer
                warning off MATLAB:timer:deleterunning
                delete(timer_obj_Camera);
            end
        catch e
            if isempty(findobj('Tag',['Fig' buffername{1}]))
                disp('Figure from Camera was closed.');
            else
                disp('An error occurred during the stream viewer update: ');
                try
                    hlp_handleerror(e);
                catch
                    disp(e.message);
                end
            end
            warning off MATLAB:timer:deleterunning
            delete(timer_obj_Camera);
        end
    end

    function on_timer_MOBIlab(varargin)
        try 
            % check if the buffer is still there
            if evalin('base',['exist(''' buffername{2} ''',''var'')'])
                
                % === update buffer contents (happens in the base workspace) ===
                
                % pull a new chunk from LSL
                [chunkdata,timestamps] = inlet{2}.pull_chunk();
                assignin('base',chunktime{2},timestamps);
                assignin('base',chunkname{2},chunkdata);
                
                % append it to the stream buffer
                evalin('base',['[' buffername{2} '.smax,' buffername{2} '.data(:,1+mod(' buffername{2} '.smax:' buffername{2} '.smax+size(' chunkname{2} ',2)-1,' buffername{2} '.pnts))] = deal(' buffername{2} '.smax + size(' chunkname{2} ',2),' chunkname{2} ');']);
                
                % get the updated stream buffer
                stream = evalin('base',buffername{2});
                
                % reformat the stream buffer to contain only the current block that should be displayed
                samples_to_get = min(stream.pnts, round(stream.srate*stream.opts.timerange));
                stream.data = stream.data(:, 1+mod(stream.smax-samples_to_get:stream.smax-1,stream.pnts));
                [stream.nbchan,stream.pnts,stream.trials] = size(stream.data);
                if ~isempty(timestamps)
                    stream.xmax = max(timestamps)-lsl_local_clock(lib);
                elseif ~isfield(stream,'xmax')
                    stream.xmax = 0;
                end
                stream.xmin = stream.xmax - (stream.pnts-1)/stream.srate;
                
                % === data post-processing for plotting ===
                
                % determine channels and samples to display
                plotchans = stream.opts.channelrange + stream.opts.pageoffset*length(stream.opts.channelrange);
                if isempty(plotchans)
                    plotchans = 1:stream.nbchan;
                else
                    plotchans = intersect(1:stream.nbchan,plotchans);
                end
                plotdata = stream.data(plotchans, round(1 : stream.srate/stream.opts.samplingrate : end));
                plottime = linspace(stream.xmin,stream.xmax,size(plotdata,2));
                
                % re-reference
                if stream.opts.reref
                    plotdata = bsxfun(@minus,plotdata,mean(plotdata)); end
                if stream.opts.standardize
                    plotdata = bsxfun(@times,plotdata,1./std(plotdata')'); end
                
                % zero-mean
				tmpdata = plotdata; tmpdata(isnan(tmpdata(:))) = 0;
                plotdata = bsxfun(@minus, plotdata, mean(tmpdata,2));
                
                % arrange for plotting
                plotoffsets = (0:size(plotdata,1)-1)*stream.opts.datascale;
                plotdata = bsxfun(@plus, plotdata', plotoffsets(2));
                
                
                % === actual drawing ===
                % smoothing
                plotdata(:,1) = smooth(plotdata(:,1),50);
                %plotdata(:,2) = smooth(plotdata(:,2),50);
                
                % draw the block contents...
                if ~isempty(plotdata)
                    if ~exist('linesMOBIlab','var') || isempty(linesMOBIlab)     
                        for k = 1 : length(opts{2}.subplot)
                            linesMOBIlab{k} = plot(axMOBIlab{k},plottime,plotdata(:,k));
                            ylabel(axMOBIlab{k},'Microvolts','FontSize',12);
                            grid(axMOBIlab{k},'minor');
                        end
                        title(axMOBIlab{1},opts{2}.streamname);
                        xlabel(axMOBIlab{length(opts{2}.subplot)},'Time (sec)','FontSize',12);
                    else
                        for k = 1 : length(linesMOBIlab)
                            set(linesMOBIlab{k},'Ydata',plotdata(:,k));
                            set(linesMOBIlab{k},'Xdata',plottime);
                        end
                    end
                
                    % update the axis limit and tickmarks
                    for k = 1 : length(linesMOBIlab)
                        axis(axMOBIlab{k},[stream.xmin stream.xmax -stream.opts.datascale size(plotdata,2)*stream.opts.datascale + stream.opts.datascale]);
                        set(axMOBIlab{k}, 'YTick',plotoffsets(2), 'YTickLabel', stream.chanlocs(k).labels);
                    end
                end
                
                drawnow;
                
            else
                try 
                    disp(['Deleting timer ' get(timer_obj_MOBIlab,'Tag') '.']);
                catch 
                    disp('Deleting timer.');
                end
                % delete the timer
                warning off MATLAB:timer:deleterunning
                delete(timer_obj_MOBIlab);
            end
        catch e
            if isempty(findobj('Tag',['Fig' buffername{2}]))
                disp('Figure from MOBIlab was closed.');
            else
                disp('An error occurred during the stream viewer update: ');
                try
                    hlp_handleerror(e);
                catch
                    disp(e.message);
                end
            end
            warning off MATLAB:timer:deleterunning
            delete(timer_obj_MOBIlab);
        end
    end

    function on_key_Camera(key)
        stream = evalin('base',buffername{1});
        switch lower(key)
            case 'uparrow'
                % decrease datascale
                stream.opts.datascale = stream.opts.datascale*0.9;
            case 'downarrow'
                % increase datascale
                stream.opts.datascale = stream.opts.datascale*1.1;
            case 'rightarrow'
                % increase timerange
                stream.opts.timerange = stream.opts.timerange*1.1;                
            case 'leftarrow'
                % decrease timerange
                stream.opts.timerange = stream.opts.timerange*0.9;                
        end
        assignin('base',buffername{1},stream);
    end

    function on_key_MOBIlab(key)
        stream = evalin('base',buffername{2});
        switch lower(key)
            case 'uparrow'
                % decrease datascale
                stream.opts.datascale = stream.opts.datascale*0.9;
            case 'downarrow'
                % increase datascale
                stream.opts.datascale = stream.opts.datascale*1.1;
            case 'rightarrow'
                % increase timerange
                stream.opts.timerange = stream.opts.timerange*1.1;                
            case 'leftarrow'
                % decrease timerange
                stream.opts.timerange = stream.opts.timerange*0.9;                
        end
        assignin('base',buffername{2},stream);
    end

    % === utility functions ===
    
    % find names of streams on the lab network...
    function names = find_streams
        streams = lsl_resolve_all(lib,0.3);
        names = cellfun(@(s)s.name(),streams ,'UniformOutput',false);
        if isempty(names)
            error('There is no stream visible on the network.'); end
    end


    % create an inlet to read from the stream with the given name
    function inlet = create_inlet(opts)
        % look for the desired device
        result = {};
        disp(['Looking for a data stream from : ' opts.streamname ' ...']);
        while isempty(result)
            result = lsl_resolve_byprop(lib,'name',opts.streamname); end
        % create a new inlet
        disp('Opening the inlet ...');
        inlet = lsl_inlet(result{1},opts.bufferrange);
    end


    % create a new stream buffer in the base workspace
    function stream = create_streambuffer(opts,info)
        stream.srate = info.nominal_srate();                                % sampling rate in Hz
        stream.chanlocs = struct('labels',derive_channel_labels(info));     % struct with per-channel meta-data
        stream.pnts = max(opts.bufferrange*stream.srate,100);               % number of data points in the buffer
        stream.nbchan = info.channel_count();                               % number of channels in the buffer
        stream.trials = 1;                                                  % number of segments in the buffer (always 1)
        stream.data = zeros(stream.nbchan,stream.pnts,stream.trials);       % the circular buffer storage
        stream.smax = 0;                                                    % number of samples that have been written into the buffer so far (wrapping around)
        stream.opts = opts;                                                 % current display options for this stream
    end


    % derive a list of channel labels for the given stream info
    function channels = derive_channel_labels(info)
        channels = deal({});
        k = 1;
        ch = info.desc().child('channels').child('channel');
        while ~ch.empty()
            name = ch.child_value_n('label');
            if name
                channels{k} = name;
                k = k + 1;
            end
            ch = ch.next_sibling_n('channel');
        end
        if length(channels) ~= info.channel_count()
            disp('The number of channels in the stream does not match the number of labeled channel records. Using numbered labels.');
            channels = cellfun(@(k)['Ch' num2str(k)],num2cell(1:info.channel_count(),1),'UniformOutput',false);
        end
    end

end

