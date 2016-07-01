function [streamsLSL] = plot_xdf()
% Plot an XDF file.
%
% In:
%   Pop-up window to choose file to plot (*.xdf or *.xdfz)
%
% Out:
%   streamsLSL : cell array of structs, one for each stream. Also a .mat 
%                file with the same name as the XDF file containing the
%                streamsLSL variable.

[filename, filepath] = uigetfile('*.xdf;*.xdfz', 'Choose an XDF file');
[streams, ~] = load_xdf(strcat(filepath,filename));

% Get entire duration of recording
for k = 1 : length(streams)
    if strcmp(streams{k}.info.name,'Length')
        Recording = k;
        RecLength = length(streams{k}.time_stamps);
        break;
    end
end

% Empty Camera
for k = 1 : length(streams)
    if strcmp(streams{k}.info.name,'Camera')
        if ~isfield(streams{k},'segments')
            streams(k) = [];
            break
        end
    end
end

% Create figures and time vectors
for k = 1 : length(streams)
    if strcmp(streams{k}.info.name,'Camera')
        streams{k}.time_seconds = streams{k}.time_stamps;
        streams{k}.time_seconds(:) = streams{k}.time_stamps(1,:) - streams{k}.time_stamps(1,1);
        Camera = k;
        axCamera = deal({});
        if (length(streams{k}.segments) > 1)
            %for i = 1 : length(streams{k}.segments)
                timeCamera = linspace(0, RecLength, (RecLength*64));
                %CameraOffset = round((streams{Camera}.segments(1).t_begin - streams{Recording}.time_stamps(1,1)),0) * 64;
                CameraOffset = round((streams{Recording}.time_stamps(1,1) - streams{Camera}.segments(1).t_begin),0) * 64;
            %end
        else
            timeCamera = linspace(0, RecLength, streams{k}.segments.num_samples);
        end
        figure('Name','Camera','NumberTitle','off');
        for i = 1 : 6
            axCamera{i} = subplot(6,1,i);
            hold(axCamera{i},'on')
            ylabel(axCamera{i},streams{k}.info.desc.channels.channel{1,i}.label,'FontSize',12);
            grid(axCamera{i},'minor');
            axCamera{i}.XLim = [0,RecLength];
        end
        title(axCamera{1},'Camera');
        xlabel(axCamera{6},'Time [sec]','FontSize',12);
    end
    if strcmp(streams{k}.info.name,'MOBIlab')
        streams{k}.time_seconds = streams{k}.time_stamps;
        streams{k}.time_seconds(:) = streams{k}.time_stamps(1,:) - streams{k}.time_stamps(1,1);        
        MOBIlab = k;
        axMOBIlab = deal({});
        if (length(streams{k}.segments) > 1)
            %for i = 1 : length(streams{k}.segments)
                timeMOBIlab = linspace(0, RecLength, (RecLength*256));
                %MOBIlabOffset = round((streams{MOBIlab}.segments(1).t_begin - streams{Recording}.time_stamps(1,1)),0) * 256;
                MOBIlabOffset = round((streams{Recording}.time_stamps(1,1) - streams{MOBIlab}.segments(1).t_begin),0) * 256;
                for kk = 1 : 2 % SMOOTHING
                    streams{MOBIlab}.time_series(kk,:) = smooth(streams{MOBIlab}.time_series(kk,:),50);
                end %
            %end
        else
            timeMOBIlab = linspace(0, RecLength, streams{k}.segments.num_samples);
            for i = 1 : 2 % SMOOTHING
                streams{MOBIlab}.time_series(i,:) = smooth(streams{MOBIlab}.time_series(i,:),50);
            end %
        end
        figure('Name','MOBIlab','NumberTitle','off');
        for i = 1 : 2
            axMOBIlab{i} = subplot(2,1,i);
            hold(axMOBIlab{i},'on')
            ylabel(axMOBIlab{i},streams{k}.info.desc.channels.channel{1,i}.label,'FontSize',12);
            grid(axMOBIlab{i},'minor');
            axMOBIlab{i}.XLim = [0,RecLength];
        end
        title(axMOBIlab{1},'MOBIlab');
        xlabel(axMOBIlab{2},'Time [sec]','FontSize',12);
    end
end

% Markers
for k = 1 : length(streams)
    if strcmp(streams{k}.info.name,'Markers')
        streams{k}.time_stamps = streams{k}.time_stamps(1,4:end);
        if exist('Camera','var')
            streams{k}.time_seconds = streams{k}.time_stamps;
            streams{k}.time_seconds(:) = streams{k}.time_stamps(1,:) - streams{Recording}.time_stamps(1,1);
            for i = 1 : 6
                minCamera = min(streams{Camera}.time_series(i,:));
                maxCamera = max(streams{Camera}.time_series(i,:));
                for ii = 1 : length(streams{k}.time_stamps)
                    line([streams{k}.time_seconds(ii) streams{k}.time_seconds(ii)],[minCamera maxCamera],'Parent',axCamera{i},'Color','r')
                end
            end
        end
        if exist('MOBIlab','var')
            streams{k}.time_seconds = streams{k}.time_stamps;
            streams{k}.time_seconds(:) = streams{k}.time_stamps(1,:) - streams{Recording}.time_stamps(1,1);
            for i = 1 : 2
                minMOBIlab = min(streams{MOBIlab}.time_series(i,:));
                maxMOBIlab = max(streams{MOBIlab}.time_series(i,:));
                for ii = 1 : length(streams{k}.time_stamps)
                    line([streams{k}.time_seconds(ii) streams{k}.time_seconds(ii)],[minMOBIlab maxMOBIlab],'Parent',axMOBIlab{i},'Color','r')
                end
            end
        end
        break
    end
end

% Draw Camera
m = 1;
if exist('Camera','var')
    if (length(streams{Camera}.segments) > 1)
        for k = 1 : 6
            for i = 1 : length(streams{Camera}.segments)
                t1 = streams{Camera}.segments(i).index_range(1) : streams{Camera}.segments(i).index_range(2);
                if t1(end) + CameraOffset > length(timeCamera)
                    CameraOffset = CameraOffset - ((t1(end) + CameraOffset) - length(timeCamera));
                end
                plot(axCamera{k},timeCamera(t1 + CameraOffset),streams{Camera}.time_series(k,t1),'b','LineWidth',0.25)
                if i+1 > length(streams{Camera}.segments)
                    break
                else
                    CameraOffset = (round((streams{Camera}.segments(i+1).t_begin - (streams{Camera}.segments(i).t_end)),0) * 64) * m;
                    m = m + 1;
                end
            end
            %CameraOffset = round((streams{Camera}.segments(1).t_begin - streams{Recording}.time_stamps(1,1)),0) * 64;
            m = 1;
            CameraOffset = round((streams{Recording}.time_stamps(1,1) - streams{Camera}.segments(1).t_begin),0) * 64;
        end
    else
        for i = 1 : 6
            plot(axCamera{i}, timeCamera, streams{Camera}.time_series(i,:))
        end
    end
end

% Draw MOBIlab
m = 1;
if exist('MOBIlab','var')
    if (length(streams{MOBIlab}.segments) > 1)
        for k = 1 : 2
            for i = 1 : length(streams{MOBIlab}.segments)
                t1 = streams{MOBIlab}.segments(i).index_range(1) : streams{MOBIlab}.segments(i).index_range(2);
                %streams{MOBIlab}.time_series(k,:) = smooth(streams{MOBIlab}.time_series(k,:),50);
                if t1(end) + MOBIlabOffset > length(timeMOBIlab)
                    MOBIlabOffset = MOBIlabOffset - ((t1(end) + MOBIlabOffset) - length(timeMOBIlab));
                end
                plot(axMOBIlab{k},timeMOBIlab(t1 + MOBIlabOffset),streams{MOBIlab}.time_series(k,t1),'Color','b','LineWidth',0.25)
                if i+1 > length(streams{MOBIlab}.segments)
                    break
                else
                    MOBIlabOffset = (round((streams{MOBIlab}.segments(i+1).t_begin - (streams{MOBIlab}.segments(i).t_end)),0) * 256) * m;
                    m = m + 1;
                end
            end
            %MOBIlabOffset = round((streams{MOBIlab}.segments(1).t_begin - streams{Recording}.time_stamps(1,1)),0) * 256;
            m = 1;
            MOBIlabOffset = round((streams{Recording}.time_stamps(1,1) - streams{MOBIlab}.segments(1).t_begin),0) * 256;
        end
    else
        for i = 1 : 2
            %streams{MOBIlab}.time_series(i,:) = smooth(streams{MOBIlab}.time_series(i,:),50);
            plot(axMOBIlab{i}, timeMOBIlab, streams{MOBIlab}.time_series(i,:))
        end
    end
end

% Create 'simplified' data
n = 1;
streamsLSL = deal({});
for k = 1 : length(streams)
    if ~strcmp(streams{k}.info.name,'Length')
        if strcmp(streams{k}.info.name,'Markers')
            streamsLSL{n}.FileName      = filename;
            streamsLSL{n}.StreamName    = streams{k}.info.name;
            streamsLSL{n}.Marks         = streams{k}.time_series;
            streamsLSL{n}.Timestamps    = streams{k}.time_seconds;
            streamsLSL{n}.AmbientLevel  = str2double(streamsLSL{n}.Marks{1}(16:end));
            streamsLSL{n}.SpeechLevel   = str2double(streamsLSL{n}.Marks{2}(15:end));
            streamsLSL{n}.TextNumber    = str2double(streamsLSL{n}.Marks{3}(14:end));
            streamsLSL{n}.Marks         = streamsLSL{n}.Marks(1,4:end);
        else
            streamsLSL{n}.StreamName    = streams{k}.info.name;
            streamsLSL{n}.Data          = streams{k}.time_series;
            streamsLSL{n}.Timestamps    = streams{k}.time_seconds;
        end
        n = n + 1;
    end
end

% varname = genvarname(filename(1:end-4));
% eval([varname ' = streamsLSL']);
% clc
% save([filename(1:end-4),'.mat'], filename(1:end-4));
