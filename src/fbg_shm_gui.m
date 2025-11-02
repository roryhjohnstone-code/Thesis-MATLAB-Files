function fbg_shm_gui()
% FBG SHM System — plots, tables, and controls (no overlap)

% ===== Figure =====
fig = figure('Name','FBG SHM System','Color','w','Position',[100 80 1200 700]);

% ===== Layout constants =====
rightX   = 0.67;         % everything on the right starts here
rightW   = 0.31;         % left plot area

% ===== Layout constants (LEFT side only) =====
leftX       = 0.06;     % left margin of the plots area
rightEdge   = 0.66;     % where plots area ends (before the right UI)
hgap        = 0.04;     % horizontal gap between columns
vgap        = 0.05;     % vertical gap between rows
topMargin   = 0.08;     % space for titles
bottomMargin= 0.08;     % space under the heatmap

colW = (rightEdge - leftX - 2*hgap)/3;  % 3 columns
rowH = 0.24;                             % height of the top/middle rows

% y origins of rows (top, middle, bottom)
row1y = 1 - topMargin - rowH;
row2y = row1y - (1.5*vgap) - rowH;
row3y = bottomMargin;                    % bottom row is the 3D heatmap

% x origins of columns
col1x = leftX;
col2x = leftX + colW + hgap;
col3x = leftX + 2*(colW + hgap);

% ===== Axes (explicit positions; prevent auto-shrink) =====
axSensors = gobjects(1,5);
axSensors(1) = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[col1x row1y colW rowH]);  title('Sensor 1 Strain');
axSensors(2) = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[col2x row1y colW rowH]);  title('Sensor 2 Strain');
axSensors(3) = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[col3x row1y colW rowH]);  title('Sensor 3 Strain');

axSensors(4) = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[col1x row2y colW rowH]);  title('Sensor 4 Strain');
axSensors(5) = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[col2x row2y colW rowH]);  title('Sensor 5 Strain');

axConsensus  = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[col3x row2y colW rowH]);  title(axConsensus,'Consensus Mean Strain & Faults');

% Bottom heatmap gets full width but shorter height so its toolbar never overlaps
ax3D = axes('Parent',fig,'Units','normalized','ActivePositionProperty','position', ...
    'Position',[leftX row3y (rightEdge-leftX) 0.20]);  % 0.20 is a good visual height
title(ax3D,'Wind Turbine Blade — Strain Heatmap');

% Consistent axes styling
for k = 1:5
    xlabel(axSensors(k),'Time (s)'); ylabel(axSensors(k),'Strain');
    grid(axSensors(k),'on'); box(axSensors(k),'on'); hold(axSensors(k),'on');
end
xlabel(axConsensus,'Time (s)'); ylabel(axConsensus,'Mean Strain');
grid(axConsensus,'on'); box(axConsensus,'on'); hold(axConsensus,'on');

% ===== Right: Fault Summary =====
uicontrol('Style','text','Parent',fig,'String','Fault Summary', ...
    'Units','normalized','Position',[rightX 0.90 rightW 0.04], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');

htFault = uitable('Parent',fig,'Units','normalized', ...
    'Position',[rightX 0.58 rightW 0.32], ...
    'ColumnName',{'Sensor','Fault Type','Start Time (s)','End Time (s)'}, ...
    'Data',cell(0,4));
setappdata(fig,'faultSummaryTableHandle', htFault);

% ===== Right: Detection Metrics =====
uicontrol('Style','text','Parent',fig,'String','Detection Metrics', ...
    'Units','normalized','Position',[rightX 0.51 rightW 0.04], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');

htMetrics = uitable('Parent',fig,'Units','normalized', ...
    'Position',[rightX 0.27 rightW 0.24], ...
    'ColumnName',{'Sensor','TP','FP','FN','Delay (s)','Accuracy (%)'}, ...
    'Data',cell(0,6));
setappdata(fig,'metricsTableHandle', htMetrics);

% After you create the Fault Summary uitable:
faultTbl.Tag = 'tblFaultSummary';
setappdata(fig,'tblFaultSummary', faultTbl);

% After you create the Detection Metrics uitable:
metricsTbl.Tag = 'tblMetrics';
setappdata(fig,'tblMetrics', metricsTbl);


% ===== Controls =====
ctrlPanel = uipanel('Parent',fig,'Title','Controls','FontWeight','bold', ...
    'Units','normalized','BackgroundColor','w','Position',[rightX 0.02 rightW 0.22]);

xL = 0.05; wL = 0.48; row = @(k) 0.86 - 0.13*(k-1);

cbFaults = uicontrol(ctrlPanel,'Style','checkbox','String','Inject Faults', ...
    'Units','normalized','Position',[xL row(1) wL 0.10], 'Value',1,'BackgroundColor','w');

cbML = uicontrol(ctrlPanel,'Style','checkbox','String','Use ML', ...
    'Units','normalized','Position',[xL row(2) wL 0.10], 'Value',1,'BackgroundColor','w');

cbKalman = uicontrol(ctrlPanel,'Style','checkbox','String','Use Kalman', ...
    'Units','normalized','Position',[xL row(3) wL 0.10], 'Value',1,'BackgroundColor','w');

cbAnim = uicontrol(ctrlPanel,'Style','checkbox','String','Animate', ...
    'Units','normalized','Position',[xL row(4) wL 0.10], 'Value',0,'BackgroundColor','w');

cbExport = uicontrol(ctrlPanel,'Style','checkbox','String','Export Animation', ...
    'Units','normalized','Position',[xL row(5) wL 0.10], 'Value',0,'BackgroundColor','w');

uicontrol(ctrlPanel,'Style','text','String','Fault type:', ...
    'Units','normalized','Position',[0.58 row(1) 0.18 0.10], ...
    'BackgroundColor','w','HorizontalAlignment','left');

popupFaultType = uicontrol(ctrlPanel,'Style','popupmenu', ...
    'Units','normalized','Position',[0.78 row(1) 0.18 0.10], ...
    'String',{'Drift','Delamination','Crack','Impact','All'}, 'Value',5);

btnRun = uicontrol(ctrlPanel,'Style','pushbutton','String','Run Simulation', ...
    'Units','normalized','Position',[0.58 0.06 0.38 0.18]);

% Run button wires into your existing runner
btnRun.Callback = @(~,~) runSimulationWithOptions( ...
    cbFaults, cbML, cbKalman, cbAnim, cbExport, popupFaultType, ...
    axSensors, axConsensus, ax3D, fig);

end


% ---------- tiny wrapper to call runSimulation ----------
function runSimulationWithOptions(cbFaults,cbML,cbKalman,cbAnim,cbExport,popupFaultType,...
                                  axSensors,axConsensus,ax3D,fig)
try
    injectFaultsFlag = logical(get(cbFaults,'Value'));
    useML            = logical(get(cbML,'Value'));
    useKalman        = logical(get(cbKalman,'Value'));
    animateFlag      = logical(get(cbAnim,'Value'));
    exportFlag       = logical(get(cbExport,'Value'));
    ftList = get(popupFaultType,'String');
    faultType = ftList{ get(popupFaultType,'Value') };

    % Call your existing simulator (unchanged signature)
    runSimulation(injectFaultsFlag, useML, useKalman, faultType, ...
                  axSensors, axConsensus, ax3D, animateFlag, exportFlag, fig);
   catch ME
    % Show full stack so we can see which line actually failed
    warning('FBG:UI:RunFailed', 'Run failed. Full report:\n%s', ...
        getReport(ME, 'extended', 'hyperlinks', 'off'));
    % Uncomment this while debugging to break at the true error site:
    % rethrow(ME)
end
end

