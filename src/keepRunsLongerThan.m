function m = keepRunsLongerThan(m, K)
    m = logical(m(:).'); if ~any(m) || K<=1, return; end
    d = diff([0 m 0]); b = find(d==1); e = find(d==-1)-1;
    keep = false(size(m));
    for i = 1:numel(b)
        if e(i)-b(i)+1 >= K, keep(b(i):e(i)) = true; end
    end
    m = keep;
end