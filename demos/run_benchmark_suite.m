function T = run_benchmark_suite()
% RUN_BENCHMARK_SUITE
% One-click grid: faults × runs × SNR × {ML,Kalman} -> CSV in /exports

%% ===== USER KNOBS =====
N_PER_FAULT = 250;                  % e.g., 10 quick, 30 solid, 50–100 thesis
SNR_LIST    = [50 30 20 15 10 5 0 -5];    % dB; Inf = no added noise
USE_ML_SET  = [true false];              % try [false true] to compare
USE_KAL_SET = [true false];              % try [false true] to compare
RNG_SEED    = 42;  % reproducible sweep
REPEATS     = 4;
FAULTS      = ["Drift","Delamination","Crack","Impact"];
rows = [];
setappdata(0,'FBG_BATCH_MODE',true);
cleanupObj = onCleanup(@() setappdata(0,'FBG_BATCH_MODE',false));

for rep = 1:REPEATS
    rng(RNG_SEED + rep);   % reproducible but different per repeat
    for f = FAULTS
        for m = USE_ML_SET
            for k = USE_KAL_SET
                for snrDB = SNR_LIST
                    cfg = struct( ...
                        'faultType',   char(f), ...
                        'useML',       logical(m), ...
                        'useKalman',   logical(k), ...
                        'snrDB',       double(snrDB), ...
                        'nPerFault',   N_PER_FAULT, ...
                        'repeat',      rep, ...
                        'seed',        RNG_SEED + rep );
                    R = simulate_once(cfg);                 % <- unchanged signature below
                    rows = [rows; pack_row(R,cfg)]; %#ok<AGROW>
                end
            end
        end
    end
end
T = struct2table(rows);
% write CSV
outDir = fullfile(pwd,'exports');
if ~exist(outDir,'dir'), mkdir(outDir); end
ts = string(datetime('now'),'yyyyMMdd_HHmmss');
outFile = fullfile(outDir, "benchmark_" + ts + ".csv");
writetable(T, outFile);
disp("✅ Benchmark CSV: " + outFile);

end % main

function R = simulate_once(cfg)
    fig = figure('Visible','off','Position',[10 10 1100 760],'Color','w');
    axS = gobjects(5,1);
    axS(1) = subplot(3,3,1,'Parent',fig);
    axS(2) = subplot(3,3,2,'Parent',fig);
    axS(3) = subplot(3,3,3,'Parent',fig);
    axS(4) = subplot(3,3,4,'Parent',fig);
    axS(5) = subplot(3,3,5,'Parent',fig);
    axC    = subplot(3,3,6,'Parent',fig);
    ax3    = subplot(3,3,[7 8 9],'Parent',fig);

    % batch hints
    setappdata(fig,'SKIP_REPORT',true);
    setappdata(fig,'SNR_DB',cfg.snrDB);
    setappdata(fig,'N_PER_FAULT',cfg.nPerFault);

    try
        runSimulation(true, cfg.useML, cfg.useKalman, cfg.faultType, ...
                      axS, axC, ax3, false, false, fig);
        R = evalin('base','simResults');
    catch ME
        warning('simulate_once:RunFailed','simulate_once failed: %s', ME.message);
        % return a minimal struct so pack_row can still proceed
        R = struct('AccuracyPerSensor',[], 'AccuracyOverall',NaN, ...
                   'ConfusionMatrix',[], 'DelaySec',[]);
    end

    if isgraphics(fig), close(fig); end
end

function row = pack_row(R,cfg)
    % ---------- config columns ----------
    row = struct();
    row.repeat     = double(cfg.repeat);
    row.seed       = double(cfg.seed);
    row.fault      = string(cfg.faultType);
    row.useML      = logical(cfg.useML);
    row.useKalman  = logical(cfg.useKalman);
    row.snrDB      = double(cfg.snrDB);
    row.nPerFault  = double(cfg.nPerFault);

    % ---------- Accuracy ----------
    accSens = safegetf(R,'AccuracyPerSensor',[]);
    accAll  = safegetf(R,'AccuracyOverall',NaN);

    if ~isempty(accSens)
        row.acc = mean(double(accSens(:)),'omitnan');   % % already, averaged across sensors
    elseif ~isnan(accAll)
        a = double(accAll);
        % If a is in [0,1], interpret as fraction and convert to %.
        % If it's already %, leave it alone.
        if a <= 1.5, a = a * 100; end
        row.acc = a;
    else
        row.acc = NaN;
    end

    % ---------- Macro precision/recall/F1 ----------
    cm = safegetf(R,'ConfusionMatrix',[]);
    if ~isempty(cm)
        cm = double(cm);
        tp = diag(cm);
        fp = sum(cm,1)' - tp;
        fn = sum(cm,2)  - tp;

        prec = mean(tp ./ max(tp+fp, eps), 'omitnan');
        rec  = mean(tp ./ max(tp+fn, eps), 'omitnan');
        f1   = 2*prec*rec / max(prec+rec, eps);
    else
        prec = NaN; rec = NaN; f1 = NaN;
    end
    row.prec_macro = prec;
    row.rec_macro  = rec;
    row.f1_macro   = f1;

    % ---------- Delay (median) ----------
    delay = safegetf(R,'DelaySec',[]);
    if ~isempty(delay)
        row.delay_med_s = median(double(delay(:)),'omitnan');
    else
        row.delay_med_s = NaN;
    end

    % ---------- Variant label ----------
    row.variant = sprintf('ML=%d,Kal=%d', row.useML, row.useKalman);
end
function v = safegetf(S, fld, def)
% SAFEGETF  robust struct getter (no arithmetic tricks)
    if isstruct(S) && isfield(S, fld) && ~isempty(S.(fld))
        v = S.(fld);
    else
        v = def;
    end
end

