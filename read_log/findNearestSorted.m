
function [idx, val, d] = findNearestSorted(cellTimes, target)
% cellTimes : sorted cell array of time strings like '15:22:36.055' or '15:22:36.0575'
% target    : time string to search
% output: idx -> nearest index, val -> nearest string, d -> absolute duration difference

    fmt = 'hh:mm:ss.SSSSSS';                          
    dt = duration(cellTimes, 'InputFormat', fmt);     
    t  = duration(target,     'InputFormat', fmt);

    n = numel(dt);
    if n == 0, error('cellTimes is empty.'); end
    if t <= dt(1), idx = 1; val = cellTimes{1}; d = dt(1) - t; return; end
    if t >= dt(n), idx = n; val = cellTimes{n}; d = t - dt(n); return; end

    lo = 1; hi = n;
    while hi - lo > 1
        mid = floor((lo + hi)/2);
        if dt(mid) < t
            lo = mid;
        else
            hi = mid;
        end
    end

    [~, pick] = min([abs(dt(lo) - t), abs(dt(hi) - t)]);
    v   = [lo hi];
    idx = v(pick);
    val = cellTimes{idx};
    d   = abs(dt(idx) - t);
end
