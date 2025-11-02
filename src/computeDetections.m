function [driftMask, delamMask, crackMask, impactMask] = computeDetections(X, t)
% COMPUTEDETECTIONS  Robust, SxN logical masks for Drift/Delam/Crack/Impact
% X: SxN signals, t: 1xN time (seconds)

    % ---------- guards ----------
    X = double(X);
    if isrow(t), t = t(:)'; end
    [S,N] = size(X);
    assert(numel(t)==N, 'computeDetections: time length must match size(X,2)');

    Fs = round((N-1) / max(t(end)-t(1), eps));
    if Fs < 2, Fs = 20; end

    % ---------- outputs ----------
    driftMask  = false(S,N);
    delamMask  = false(S,N);
    crackMask  = false(S,N);
    impactMask = false(S,N);

    % ---------- per-sensor detection ----------
    for s = 1:S
        x = X(s,:);

        % z-normalize for robust thresholds
        mu = mean(x); sg = std(x); if sg==0, sg=1; end
        xz = (x - mu) / sg;

        % ===== Impact: short, strong spike (<= 3 samples) =====
        wE   = max(3, round(0.10*Fs));               % ~0.1 s
        E    = movsum(xz.^2, wE, 'Endpoints','shrink');
        Emed = median(E); sE = 1.4826*mad(E,1);
        thrI = Emed + 4.5*sE;                           % quite strict
        mkI  = E > thrI;

        % keep only the strongest very short burst (<= 3 samples)
        mkI  = keepShortestStrongRun(mkI, 5, xz);

        % initial halo (used below by crack); will refresh after delam
        mkIwide = dilateMask(mkI, round(0.6*Fs));

        % --- CRACK (step with sustained offset, not oscillatory) ---
        minStep = 1.0;                        % ↑ larger step
        holdWin = max(round(0.60*Fs), 12);    % ↑ longer sustain
        halfWin = max(round(0.50*Fs), 10);    % ↑ wider context

        mkC = false(1,N);
        dx  = [0 diff(x)];
        [~, idx] = max(movmean(max(dx,0), round(0.15*Fs)));   % candidate up-step

        i1 = max(1, idx-halfWin);  i2 = min(N, idx+halfWin);
        pre  = i1:idx-1;  post = idx:min(N, idx+holdWin-1);

        if numel(pre) > 5 && numel(post) > 5
            stepAmp = mean(x(post)) - mean(x(pre));
            % Require: large step, and post region is not too "noisy/oscillatory"
            if stepAmp >= minStep && std(x(post)) < 0.8*std(x(pre)) + 0.10
                mkC(post) = true;
            end
        end

        % Avoid calling crack inside strong impacts
        mkC = mkC & ~mkIwide;
        mkCwide = dilateMask(mkC, round(0.6*Fs));
        crackMask(s,:) = mkC;

      % ===== Delamination: band energy burst with envelope =====
nyq = Fs/2;
% safe, mid-low band; slightly wider
f1  = max(1.0, min(0.15*nyq, 0.03*Fs));
f2  = min(0.55*nyq, max(f1+1.5, 0.12*Fs));

bp   = safeBandpass(x, Fs, [f1 f2]);
env  = movmean(abs(hilbert(bp)), max(5, round(0.20*Fs)));

muE  = median(env); sdE = 1.4826*mad(env,1);
thrD = muE + 2.4*sdE;                            % WAS 9*sdE

mkD  = env > thrD;

% require ~0.6 s hold (WAS 1.2 s)
minHoldD = max(round(0.50*Fs), 8);
mkD  = keepRunsLongerThan(mkD, minHoldD);
mkD  = keepLongestRun(mkD);
% Require oscillatory behavior inside the delam segment (not a monotonic ramp)
if any(mkD)
    dxb     = abs(diff(x));                     % baseline diff stats
    baseD   = median(dxb);
    baseMAD = 1.4826*mad(dxb,1);

    segIdx  = find(mkD);
    if numel(segIdx) >= 3
        oscScore = mean(abs(diff(x(segIdx))));
        if ~(oscScore > baseD + 2.5*baseMAD)
            mkD(:) = false;   % kill delam if segment isn't oscillatory
        end
    end
end
if any(mkD)
    segIdx  = mkD;
    envSeg  = env(segIdx);
    [pks, ~] = findpeaks(envSeg);
    if numel(pks) < 3
        mkD(:) = false;   % not "bursty" enough → likely a step/edge, not delam
    end
end
if any(mkD)
    segIdx = find(mkD);
    if numel(segIdx) >= 12
        xv   = x(segIdx);
        n    = numel(xv);
        a    = mean(xv(1:round(n/3)));
        c    = mean(xv(round(2*n/3):end));
        if (c - a) > 0.8   % step magnitude gate; tune if needed (0.6–1.2)
            mkD(:) = false;
        end
    end
end

% Fallback: local variance burst if envelope found nothing
if ~any(mkD)
    wv  = max(round(0.40*Fs), 8);
    lv  = movvar(x, wv, 0, 'Endpoints','shrink');
    muV = median(lv); sdV = 1.4826*mad(lv,1);
    mkV = lv > (muV + 2.5*sdV);
    mkV = keepRunsLongerThan(mkV, minHoldD);
    mkV = keepLongestRun(mkV);
    if any(mkV), mkD = mkV; end
end

% avoid overlaps with higher priority classes
mkD = mkD & ~mkIwide & ~mkCwide;
delamMask(s,:) = mkD;

% ===== (UPDATE) Impact suppression near delamination, then refresh halo =====
mkI  = mkI & ~dilateMask(mkD, round(0.40*Fs));   % NEW: kill impact inside/near delam
mkIwide = dilateMask(mkI, round(0.6*Fs));        % refresh halo after suppression
impactMask(s,:) = mkI;                           % store FINAL impact mask

% ===== Drift: persistent positive trend (robust + global gates) =====
winSec = 2.0;
W = max(round(winSec*Fs), 8);

slopePS = zeros(1,N);
posFrac = zeros(1,N);

for i = 1:N
    j1 = max(1, i-W+1);
    j2 = i;
    tt = t(j1:j2) - t(j1);
    yy = x(j1:j2);
    if numel(tt) >= 8
        p = polyfit(tt, yy, 1);          % slope per second
        slopePS(i) = p(1);
        posFrac(i) = mean(diff(yy) >= 0);
    end
end

% --- global trend / rise gates (quiet in non-drift tests)
pg   = polyfit(t - t(1), x, 1);          % global slope per second
rise = x(end) - x(1);                    % net rise over record
allowDrift = (pg(1) > 0.02) & (rise > max(0.12, 0.15*std(x)));

% --- local consistency + threshold
thrSlope = max(0.020, median(slopePS) + 0.6 * 1.4826 * mad(slopePS,1));
mkR = (slopePS > thrSlope);  %(posFrac > 0.60);

% persistence (be a bit stricter)
mkR = keepRunsLongerThan(mkR, max(round(0.8*Fs), 12));

% stronger veto halos around higher-priority faults
mkR = mkR & ~dilateMask(mkI, round(0.80*Fs)) ...
          & ~dilateMask(mkC, round(0.60*Fs)) ...
          & ~dilateMask(mkD, round(0.50*Fs));

% apply global gates last
if ~allowDrift
    mkR(:) = false;
end

% fallback: if absolutely nothing but clear global ramp, mark a short tail
if ~any(mkR) && allowDrift && (pg(1) > 0.04)
    idx = round(0.7*N):N;
    mkR(idx) = true;
end

driftMask(s,:) = mkR;

    end

    % ===== Priority pass (Impact > Crack > Delam > Drift) =====
    for s = 1:S
        imp = impactMask(s,:);
        crk = ~imp & crackMask(s,:);
        dlm = ~imp & ~crk & delamMask(s,:);
        drf = ~imp & ~crk & ~dlm & driftMask(s,:);
        impactMask(s,:) = imp;
        crackMask(s,:)  = crk;
        delamMask(s,:)  = dlm;
        driftMask(s,:)  = drf;
    end
end

% ----------------- helpers -----------------
function y = safeBandpass(x, Fs, band)
    try
        y = bandpass(x, band, Fs, "Steepness",0.85, "StopbandAttenuation",60);
    catch
        hp = highpass(x, max(0.2, band(1)*0.8), Fs);
        y  = lowpass(hp, min(Fs/2-1e-3, band(2)*1.1), Fs);
    end
end

function m = keepShortestStrongRun(m, maxLen, xref)
    m = logical(m(:).'); if ~any(m), return; end
    d = diff([0 m 0]); b = find(d==1); e = find(d==-1)-1;
    best = -inf; keep = false(size(m));
    for i = 1:numel(b)
        seg = b(i):e(i);
        if numel(seg) <= maxLen
            pk = max(abs(xref(seg)));
            if pk > best, best = pk; keep = false(size(m)); keep(seg) = true; end
        end
    end
    m = keep;
end

function mkOut = dilateMask(mk, w)
% 1D morphological dilation by a half-width 'w' (in samples)
    w = max(0, round(w));
    if w == 0, mkOut = mk(:).'; return; end
    k = ones(1, 2*w+1);
    mkOut = conv(double(mk(:).'), k, 'same') > 0;
end
function mOut = keepRunsPassing(mIn, okFn)
% Keep only those 1-runs for which okFn(segmentIndices) returns true.
    m = logical(mIn(:).'); 
    if ~any(m), mOut = m; return; end
    d = diff([0 m 0]); b = find(d==1); e = find(d==-1)-1;
    keep = false(size(m));
    for i = 1:numel(b)
        seg = b(i):e(i);
        if okFn(seg), keep(seg) = true; end
    end
    mOut = keep;
end
