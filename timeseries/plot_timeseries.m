function [h] = plot_timeseries(timeseries, varargin)
% Plot time series. timeseries is a struct array containing three fields:
% timeseries.var: Variable name
% timeseries.times: An nx1 array of Matlab datetimes reflecting the
% timestamps
% timeseries.vals: An nx1 array of values corresponding to the timestamps
% in timeseries.times
% Can take in optional cell array of strings to search for so that we only
% plot those.
% Can additionally take in optional cell array of strings to search for
% when plotting text (e.g. med admin values).

if(nargin > 1)
    filter = varargin{1};
else
    filter = {};
end
if(nargin > 2)
    text_filter = varargin{2};
else
    text_filter = {};
end

h = figure;
hold on;

legendvars = {};
for i = 1:length(timeseries)
    if(~isempty(timeseries(i).times))
        if(isempty(filter) || contains(timeseries(i).var,filter,'IgnoreCase',true))
            [times_sorted, sort_inds] = sort(timeseries(i).times);
            vals_sorted = timeseries(i).vals(sort_inds);
            if(isnumeric(vals_sorted))
                plot(times_sorted, vals_sorted, '.-','DisplayName',sprintf('(%d) %s',i,timeseries(i).var));
            elseif(contains(timeseries(i).var,'StartStop'))
                for j = 1:2:length(times_sorted)-1
                    plot(times_sorted(j:j+1),[1 1],'+-','LineWidth',2,'DisplayName',timeseries(i).var);
                end
            elseif(iscellstr(vals_sorted))
                y_counter = 0;
                for j = 1:length(times_sorted)
                    if(isempty(text_filter) || contains(vals_sorted{j}, text_filter,'IgnoreCase',true))
                        y_counter = y_counter + 1;
                        th = text(times_sorted(j),50*mod(y_counter,2),vals_sorted{j},'Rotation',90,'FontSize',6);
                        uistack(th,'bottom');
                    end
                end
            end
        end
    end
end
legend();
yl = ylim;
ylim([0, yl(2)]);