function m = keepLongestRun(m)
    m = logical(m(:).'); if ~any(m), return; end
    d = diff([0 m 0]); b = find(d==1); e = find(d==-1)-1;
    [~,k] = max(e-b+1);
    keep = false(size(m)); keep(b(k):e(k)) = true;
    m = keep;
end
