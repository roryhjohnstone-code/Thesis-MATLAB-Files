function plotSensorStrain(ax, t, y, driftM, ~, delamM, crackM, impactM, sIdx)
if nargin < 9 || isempty(sIdx), sIdx = 1; end
if isempty(ax) || ~isgraphics(ax), return; end
axes(ax); cla(ax); hold(ax,'on'); grid(ax,'on');

N = numel(t);
pad = @(m) logical(reshape(local_padMask(m, N), 1, N));
driftM   = pad(driftM);
delamM   = pad(delamM);
crackM   = pad(crackM);
impactM  = pad(impactM);

% base strain
plot(t, y, 'b-', 'LineWidth', 1.2, 'DisplayName','Strain');

% markers (no legend entries for per-sensor, keeps it clean)
if any(driftM)
    plot(t(driftM),  y(driftM),  'r.',  'MarkerSize', 12, 'HandleVisibility','off');
end
if any(delamM)
    plot(t(delamM), y(delamM), 'gs', 'MarkerSize', 5, 'MarkerFaceColor',[0.3 0.8 0.3], ...
         'LineStyle','none', 'HandleVisibility','off');
end
if any(crackM)
    d  = diff([false crackM false]);
    b  = find(d==1, 1, 'first');
    if ~isempty(b)
        plot(t(b), y(b), 'kd', 'MarkerSize', 7, 'MarkerFaceColor','k', 'HandleVisibility','off');
    end
end
if any(impactM)
    plot(t(impactM), y(impactM), 'p', 'MarkerSize', 9, 'LineStyle','none', ...
         'Color',[0 0.6 0.6], 'MarkerFaceColor',[0 0.6 0.6], 'HandleVisibility','off');
end

xlabel('Time (s)'); ylabel('Strain');
title(sprintf('Sensor %d Strain', sIdx));
set(ax,'Box','off');   % <- removes the black edge box look
end


