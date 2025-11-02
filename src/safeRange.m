function idx = safeRange(a,b,N)
%SAFERANGE  Build a:b safely, clamped to [1..N] and integer.
% Returns [] if inputs aren't valid scalars.
    if nargin < 3 || isempty(N), N = inf; end
    a = round(a); b = round(b);
    if ~isscalar(a) || ~isscalar(b) || ~isfinite(a) || ~isfinite(b)
        idx = []; return;
    end
    a = max(1, min(N, a));
    b = max(1, min(N, b));
    if a > b, [a,b] = deal(b,a); end
    idx = a:b;
end
