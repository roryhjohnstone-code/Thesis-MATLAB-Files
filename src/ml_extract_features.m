function f = ml_extract_features(x,Fs)
% Return row vector of features from 1-D signal x
x = x(:)'; N = numel(x);
% Detrend gently
xd = detrend(x,1);

% --- basic stats
m  = mean(xd); s = std(xd); r = range(xd);
sk = skewness(xd); ku = kurtosis(xd);
zcr = sum(abs(diff(sign(xd))))/N;

% --- slope-ish (drift cue)
p = polyfit(1:N, xd, 1); slope = p(1);

% --- Welch spectrum (robust to very short windows)
L = numel(xd);
if L < 8
    % If the slice is tiny, skip spectral stats gracefully
    Px = 1; F = 0;        % dummy spectrum
else
    win   = min(L, max(16, round(0.8*Fs)));   % window ≤ signal length
    nover = max(0, min(round(0.5*win), win-1));
    [Px,F] = pwelch(xd, win, nover, [], Fs, 'psd');
    Px = Px(:)'; F = F(:)';
    Px = Px / (sum(Px)+eps);
end


% spectral centroid & spread
sc = sum(F.*Px);
ss = sqrt(sum(((F - sc).^2).*Px));

% band energies (0–0.6, 0.6–1.5, 1.5–3 Hz)
b1 = bandpower(Px,F,[0 0.6],'psd');
b2 = bandpower(Px,F,[0.6 1.5],'psd');
b3 = bandpower(Px,F,[1.5 3.0],'psd');

% --- short-time energy for impacts
w = round(0.15*Fs);
E = movsum(xd.^2, w, 'Endpoints','shrink');
Epk = max(E);  Emed = median(E);
impulseRatio = Epk / (Emed + eps);

% --- step detector cue (crack)
dx = diff(xd);
stepStat = max(0, max(movmean(dx,round(0.15*Fs))) - median(dx));

% --- delam cue: local variance burst
locVar = movvar(xd, round(0.4*Fs), 0, 'Endpoints','shrink');
delamStat = prctile(locVar,95);

f = [m s r sk ku zcr slope sc ss b1 b2 b3 impulseRatio stepStat delamStat];
end
