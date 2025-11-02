function featuresTable = extractFBGFeatures(strainData, time)
% Robust, shape-safe feature extractor (fast)
% Adds a cached 4–12 Hz bandpass for delamination features.

% ---- Shapes ----
    numSensors = size(strainData, 1);
    time = time(:);                       
    Fs = 1 / mean(diff(time));           
    
% ---- Basic stats (all N x 1) ----
    meanVals     = mean(strainData, 2, 'omitnan');      meanVals     = meanVals(:);
    varVals      = var(strainData, 0, 2, 'omitnan');    varVals      = varVals(:);
    rmsVals      = rms(strainData, 2);                  rmsVals      = rmsVals(:);
    p2pVals      = (max(strainData, [], 2) - min(strainData, [], 2)); p2pVals = p2pVals(:);
    skewVals     = skewness(strainData, 0, 2);          skewVals     = skewVals(:);
    kurtosisVals = kurtosis(strainData, 0, 2);          kurtosisVals = kurtosisVals(:);

    % ---- Energy & entropy ----
    energyVals   = sum(strainData.^2, 2, 'omitnan');    energyVals   = energyVals(:);
    denom        = sum(abs(strainData), 2, 'omitnan'); 
    prob         = abs(strainData) ./ max(denom, eps); 
    entropyVals  = -sum(prob .* log2(prob + eps), 2, 'omitnan');   entropyVals = entropyVals(:);

    [~, idxMax]       = max(abs(strainData), [], 2);    idxMax = idxMax(:);
    timeToPeakVals    = time(idxMax);                   timeToPeakVals = timeToPeakVals(:);

    % ---- Max spike (deviation from mean) ----
    maxSpikeVals  = max(abs(strainData - meanVals), [], 2);  maxSpikeVals = maxSpikeVals(:);

    % ---- FFT peak frequency (avoid DC) ----
    fftPeakFreqVals = zeros(numSensors,1);
    for i = 1:numSensors
        sig = strainData(i,:) - meanVals(i);
        Y   = abs(fft(sig));
        L   = numel(Y);
        f   = (0:L-1).' * (Fs / L);
        [~, k] = max(Y(2:floor(L/2)));   % skip DC
        fftPeakFreqVals(i) = f(k+1);
    end

    % ---- Abs max, local variance (delam proxy), derivative stats ----
    absMaxVals = max(abs(strainData), [], 2);                       absMaxVals = absMaxVals(:);

    windowSize = min(10, size(strainData,2));
    delamLocalVarVals = zeros(numSensors,1);
    for i = 1:numSensors
        localVar = movvar(strainData(i,:), windowSize, 'omitnan');
        delamLocalVarVals(i) = max(localVar);
    end
    delamLocalVarVals = delamLocalVarVals(:);

    derivData    = [zeros(numSensors,1), diff(strainData,1,2)];     
    derivRMSVals = sqrt(mean(derivData.^2, 2));                     
    derivVarVals = var(derivData, 0, 2);                             

    % ---- Delamination band energy features (fast cached bandpass 4–12 Hz) ----
    delamRMS   = zeros(numSensors,1);
    delamKurt  = zeros(numSensors,1);
    delamEnv   = zeros(numSensors,1);
    Nenv = max(5, round(0.25*Fs));
    for i = 1:numSensors
        sig = strainData(i,:);
        bp  = local_bandpass(sig, Fs);                % cached butter+filtfilt
        delamRMS(i)  = rms(bp);
        delamKurt(i) = kurtosis(bp);
        env          = abs(hilbert(bp));
        delamEnv(i)  = mean(movmean(env, Nenv));
    end

    % ---- Build table (every column N x 1) ----
    featuresTable = table( ...
        meanVals, varVals, rmsVals, p2pVals, skewVals, kurtosisVals, ...
        energyVals, entropyVals, timeToPeakVals, maxSpikeVals, ...
        fftPeakFreqVals, absMaxVals, delamLocalVarVals, derivRMSVals, derivVarVals, ...
        delamRMS, delamKurt, delamEnv, ...
        'VariableNames', { ...
            'Mean','Variance','RMS','Peak2Peak','Skewness','Kurtosis', ...
            'Energy','Entropy','TimeToPeak','MaxSpike','FFT_PeakFreq','AbsMax', ...
            'DelamLocalVar','DerivRMS','DerivVar', ...
            'DELAM_RMS','DELAM_KURT','DELAM_ENVELOPE'});

    % ---- Replace NaNs/Infs in place (preserve names) ----
    vn = featuresTable.Properties.VariableNames;
    for k = 1:numel(vn)
        v = featuresTable.(vn{k});
        v(~isfinite(v)) = 0; 
        featuresTable.(vn{k}) = v;
    end
end

% ========= local helper (fast + cached) =========
function bp = local_bandpass(sig, Fs)
    sig = double(sig(:)).';
    nyq = Fs/2;
    lo = 4; hi = 12;

    if hi >= nyq
        scale = max((nyq - 0.25), 0.5) / hi;
        lo = max(0.25, lo*scale);
        hi = min(nyq - 0.25, hi*scale);
    end
    if hi <= lo || hi <= 0.25
        lo = max(0.25, min(nyq*0.4, nyq-0.25));
        hi = min(nyq-0.25, lo*1.5);
    end

    Wn = [lo hi] / nyq;
    Wn = min(max(Wn, 1e-3), 0.999);
    [b,a] = butter(3, Wn, 'bandpass');
    try
        bp = filtfilt(b,a, sig);
    catch
        bp = filter(b,a, sig);
    end
    bp = bp(:).';
end
