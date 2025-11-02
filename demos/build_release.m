function build_release()
% Build a clean /release folder with only necessary files for examiners.
% - Copies deps for your entry points into release/src (both .m and .mat)
% - Always copies demos/*.m into release/demos (including run_benchmark_suite.m)
% - Writes MATLAB_TOOLBOXES.md, README.md, .gitignore

%% --- Locate paths ---
toolsDir   = fileparts(mfilename('fullpath'));
projectDir = fileparts(toolsDir);
relDir     = fullfile(projectDir, 'release');
srcOut     = fullfile(relDir, 'src');
demoOut    = fullfile(relDir, 'demos');
dataOut    = fullfile(relDir, 'data'); %#ok<NASGU>
docsOut    = fullfile(relDir, 'docs'); %#ok<NASGU>

%% --- Entry points (adjust if needed) ---
entryPoints = { ...
    'fbg_shm_gui.m', ...
    'runSimulation_v2.m' ...
    % (We will copy demos explicitly below; no need to add them here)
};

%% --- Put project on path and get dependencies ---
addpath(genpath(projectDir));
deps = listDependencies(entryPoints);

%% --- Clean /release and recreate folders ---
if exist(relDir,'dir'), rmpath(genpath(relDir)); end   % avoid "Removed … from path" spam
if exist(relDir,'dir'), rmdir(relDir,'s'); end
mkdir(srcOut); mkdir(demoOut); mkdir(fullfile(relDir,'docs')); mkdir(fullfile(relDir,'data'));

%% --- Copy only your project files (avoid MATLAB internals) ---
isLocal    = contains(deps.files, projectDir);
keepExt    = endsWith(deps.files, {'.m','.mat'}, 'IgnoreCase', true);
filesToCopy = deps.files(isLocal & keepExt);

copied = 0;
for i = 1:numel(filesToCopy)
    [~, name, ext] = fileparts(filesToCopy{i});
    dest = fullfile(srcOut, [name ext]);
    copyfile(filesToCopy{i}, dest);
    copied = copied + 1;
end

%% --- Write toolbox list ---
fid = fopen(fullfile(relDir,'MATLAB_TOOLBOXES.md'),'w'); assert(fid>0);
fprintf(fid,"# Required MATLAB Toolboxes\n\n");
for i=1:numel(deps.products)
    fprintf(fid,"- %s (version %s)\n", deps.products(i).Name, deps.products(i).Version);
end
fclose(fid);

%% --- Minimal demo runner (always present) ---
fid = fopen(fullfile(demoOut,'run_demo.m'),'w'); assert(fid>0);
fprintf(fid, "%% Minimal demo: launches GUI if present, else runs headless simulation\n");
fprintf(fid, "rng(1);\n");
fprintf(fid, "addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','src')));\n");
fprintf(fid, "if exist('fbg_shm_gui.m','file'), fbg_shm_gui; else, runSimulation_v2; end\n");
fclose(fid);

%% --- COPY ALL demos/*.m explicitly (this is where run_benchmark_suite.m is copied) ---
demoSrc = fullfile(projectDir, 'demos');
copiedDemos = {};
if exist(demoSrc, 'dir')
    d = dir(fullfile(demoSrc, '*.m'));
    for i = 1:numel(d)
        srcFile = fullfile(demoSrc, d(i).name);
        dstFile = fullfile(demoOut, d(i).name);
        copyfile(srcFile, dstFile);
        copiedDemos{end+1} = d(i).name; %#ok<AGROW>
    end
else
    fprintf('[build_release] No demos/ folder found at %s\n', demoSrc);
end

%% --- README ---
fid = fopen(fullfile(relDir,'README.md'),'w'); assert(fid>0);
fprintf(fid, "# FBG-Based SHM (MATLAB)\n\n");
fprintf(fid, "## How to run\n");
fprintf(fid, "```matlab\ncd demos\nrun_demo\n```\n\n");
fprintf(fid, "### Benchmark suite\n");
fprintf(fid, "If present:\n```matlab\ncd demos\nrun_benchmark_suite\n```\n");
fprintf(fid, "CSV will be written under `exports/`.\n");
fclose(fid);

%% --- .gitignore ---
fid = fopen(fullfile(relDir,'.gitignore'),'w'); assert(fid>0);
fprintf(fid, "## MATLAB\n*.asv\n*.m~\n*.mlx.autosave\n*.mex*\n*.mat\n*.fig\n*.slxc\n.DS_Store\n");
fclose(fid);

%% --- Report ---
fprintf('\n✅ Built release with %d source files at %s\n', copied, relDir);
if ~isempty(copiedDemos)
    fprintf('   Demos copied (%d): %s\n', numel(copiedDemos), strjoin(copiedDemos, ', '));
    % Assert the benchmark file exists if you expect it:
    if any(strcmpi(copiedDemos,'run_benchmark_suite.m'))
        fprintf('   ✓ Found and copied run_benchmark_suite.m\n');
    else
        fprintf('   ⚠ run_benchmark_suite.m not found in demos/ (skipped)\n');
    end
else
    fprintf('   ⚠ No demo .m files were copied.\n');
end
% Final step: force-copy demos (handles demos/ or demo/)
if exist(fullfile(projectDir,'tools','copy_demos.m'),'file')==2
    copy_demos(projectDir, relDir);
end

end
