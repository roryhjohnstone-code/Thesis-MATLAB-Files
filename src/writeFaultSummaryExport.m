function writeFaultSummaryExport(fig, outDir, faultType)
% Export Fault Summary uitable data to CSV (with mismatch column)

if nargin<3, faultType=""; end
if nargin<2 || isempty(outDir), outDir = fullfile(pwd,'exports'); end
if ~exist(outDir,'dir'), mkdir(outDir); end

ht = getappdata(fig,'tblFaultSummary');
if isempty(ht) || ~isvalid(ht), return; end
data = ht.Data;
if isempty(data), return; end

% Ensure all rows are cellstr
for c = 1:size(data,2)
    if isnumeric(data{1,c})
        data(:,c) = cellfun(@num2str, data(:,c), 'UniformOutput', false);
    end
end

% Add mismatch column if ML Prediction exists
hasML = size(data,2) >= 5;
if hasML
    mismatch = repmat("No", size(data,1), 1);
    for r = 1:size(data,1)
        det = string(data{r,2});
        ml  = string(data{r,5});
        if ml~="" && ml~="Normal" && ml~=det && det~="None"
            mismatch(r) = "Yes";
        end
    end
    data(:,6) = cellstr(mismatch);
    colNames = {'Sensor','Detector Fault','Start Time','End Time','ML Prediction','Mismatch'};
else
    colNames = {'Sensor','Detector Fault','Start Time','End Time'};
end

% Build table
T = cell2table(data, 'VariableNames', colNames);

% Filename with timestamp + fault type
ts = string(datetime('now'),'yyyyMMdd_HHmmss');
fname = "faultSummary_" + faultType + "_" + ts + ".csv";

writetable(T, fullfile(outDir, fname));
fprintf('✅ Fault summary exported: %s\n', fullfile(outDir, fname));
end
