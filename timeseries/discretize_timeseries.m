function [reg_vals, reg_times, vars] = discretize_timeseries(timeseries, timestep)
% Discretize time series data into uniform discrete timesteps.  timeseries
% is a struct array containing three fields:
% timeseries.var: Variable name
% timeseries.times: An nx1 array of Matlab datetimes reflecting the
% timestamps
% timeseries.vals: An nx1 array of values corresponding to the timestamps
% in timeseries.times

times = [];
for i = 1:length(timeseries)
    times = [times; timeseries(i).times(:)];
end
times = sort(times);

start_time = dateshift(min(times),'start','hour');
end_time = dateshift(max(times),'start','hour');

reg_times = start_time:timestep:end_time;
first_time = sum(times(1) >= reg_times);
reg_times = reg_times(first_time:end);
reg_vals = nan(length(timeseries),length(reg_times));

vars = {timeseries.var};

for i = 1:length(timeseries)
    if(isempty(timeseries(i).vals))
        continue;
    end
    for j = 1:length(reg_times)
        mask = timeseries(i).times >= reg_times(j) & timeseries(i).times < reg_times(j) + timestep;
        ts = timeseries(i).times(mask);
        vs = timeseries(i).vals(mask);
        if(isempty(vs))
            continue;
        elseif(length(vs) == 1)
            vint = vs;
        else
            vint = mean(vs);
        end
        reg_vals(i,j) = vint;
    end
end