function [ts] = get_ts(pt_data, ts_name)
ts = pt_data.timeseries(contains({pt_data.timeseries.var},ts_name));