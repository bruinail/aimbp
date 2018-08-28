function [variances] = calc_segmented_variance(timeseries)
    [vs, ts, ~] = discretize_timeseries(timeseries,hours(0.25));
    vs = vs';
    ts = ts';
    t_hours = hours(ts - ts(1));
    seg_length = 6;
    variances = [];
    for segs = 1:ceil(t_hours(end)/seg_length)
        valid = t_hours >= (segs-1)*seg_length & t_hours < segs*seg_length & ~isnan(vs);
        if(sum(valid) > 5)
            v_seg = vs(valid);
            t_seg = t_hours(valid);
            fit_seg = [ones(size(t_seg)) t_seg] \ v_seg;
            fit_v = fit_seg(2)*t_seg + fit_seg(1);
            variances = [variances nanvar(v_seg - fit_v)];
        end
    end
end