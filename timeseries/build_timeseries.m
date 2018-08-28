function [timeseries] = build_timeseries(chartdata,itemids,varnames)
% Builds a timeseries structure array of the following format from chart
% data and a list of itemids (corresponding to varnames):
% timeseries.var: Variable name
% timeseries.times: An nx1 array of Matlab datetimes reflecting the
% timestamps
% timeseries.vals: An nx1 array of values corresponding to the timestamps
% in timeseries.times

dt_format = 'yyyy-MM-dd HH:mm:ss.S';

for i = 1:length(itemids)
    mask = chartdata.itemid == itemids(i);
    charttimes = datetime(chartdata.charttime(mask),'InputFormat',dt_format);
    vals = chartdata.valuenum(mask);
    timeseries(i) = struct('var',varnames{i},'times',charttimes,'vals',vals);
end