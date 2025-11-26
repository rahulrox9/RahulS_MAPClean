%% MAPClean – Microstructurally Adaptive Pixel-Level Cleaning)
%
clc; clear; close all;
addpath(genpath(pwd)); 
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

%% === START OF USER CONFIGURATION === %%

%% --- 1A. Stage Control (Flags) ---
% Set flags to control which cleaning steps are executed or loaded from checkpoints
runStart    = true;    % Initial plots
runMAD      = true;    % Mean Angular Deviation Filter
runCrop     = true;    % Sample Mask/Cropping
runPhaseWSR = true;    % Phase Wild Spike Removal
runOriWSR   = true;    % Orientation Wild Spike Removal
runHoleFill = true;    % Standard Hole Filling (BFS/MPF)
runProFill  = true;    % Protected Pixel Filling (Protected Holes)
runSaveFile = true;    % Export final EBSD and parameters

%% --- 1B. Parameters (Global Configuration) ---
% NOTE: The 'global' keyword is used to access this struct in all functions.
global params
disp('Initialising Parameters...');

params.exportRes      = 300;        % dpi for figure export resolution
params.madThreshold   = 0.9;        % radians, MAD filter threshold

% Phase WSR parameters
params.radius_phase   = 2;          % WSR kernel radius (in pixels)
params.min_neighbours = 3;          % Minimum number of neighbours (unused in doPhaseWSR logic, kept for consistency)
params.min_dom_frac   = 0.5;        % Minimum dominant phase fraction for 'aggressive' flip
params.numPasses      = 25;         % Number of passes (unused, kept for consistency)

% Orientation WSR parameters
params.misTol_ori     = 5*degree;   % Orientation misorientation tolerance for clustering
params.thresholdFrac  = 0.75;       % Minimum dominant cluster fraction for WSR (Strict)
params.minLead        = 2;          % Minimum lead count for WSR (Relaxed/Aggressive)
params.scaleLead      = 0.1;        % Scaling factor for required lead (Relaxed/Aggressive)
params.minFrac_ori    = 0.25;       % Minimum fraction of similar neighbours to avoid removal (Pre-filter)
params.radius_ori     = 2;          % WSR kernel radius (in pixels)

% Hole-filling parameters
params.radius_fill    = [6 5 4 3 2 1]; % Radii used for sequential filling
% Map(radius) -> [Ni_thresh, fracDom_thresh] 
params.phaseFrac    = containers.Map('KeyType','double','ValueType','any');
params.phaseFrac(6) = [0.4 0.75];
params.phaseFrac(5) = [0.4 0.75];
params.phaseFrac(4) = [0.4 0.75];
params.phaseFrac(3) = [0.4 0.75];
params.phaseFrac(2) = [0.4 0.75];
params.phaseFrac(1) = [0.4 0.75];

% Protected Hole-filling parameters
params.thresholdFracRing  = 2/3;    % Minimum dominant cluster fraction in a ring
params.coverageFrac       = 2/3;    % Minimum indexed neighbour coverage in the kernel

disp('✔ Parameters initialised');

%% === END OF USER CONFIGURATION === %%

%% --- 2. Directories and Setup ---
dataDir       = fullfile(pwd,'DataFiles');
fileList      = dir(fullfile(dataDir, '*.ctf'));
fileList      = fileList(~contains({fileList.name}, '_')); 
checkpointDir = fullfile(pwd,'checkpoints');
exportDir     = fullfile(pwd,'exports');
% Create directories if they do not exist
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end
if ~exist(exportDir,'dir'), mkdir(exportDir); end

%% --- 3. Loop over Samples ---
% Iterate through all specified CTF files
for fi = 1:numel(fileList)
    % Extract sample name and set up output paths
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    exportPath = fullfile(exportDir, sampleName);
    if ~exist(exportPath,'dir'), mkdir(exportPath); end
    
    % Setup logging/diary file
    diaryFile = fullfile(exportDir, [sampleName '_logfile.txt']);
    if exist(diaryFile, 'file')
        delete(diaryFile); % Delete the file if it already exists
    end
    diary(diaryFile); % Start a fresh diary
    diary on;
    fprintf('\n===== Processing Sample: %s =====\n', sampleName);
    
    %% --- 3.1. Load EBSD Data ---
    ebsd_raw = EBSD.load(fullfile(dataDir,fileList(fi).name), ...
                         'convertSpatial2EulerReferenceFrame').gridify;
    % Identify phases
    phases = ebsd_raw.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    
    %% --- 3.2. Set Phase Colours ---
    % Custom colour assignment for specific minerals
    ForsteriteColor = str2rgb('red');
    DiopsideColor   = str2rgb('blue');
    AnorthiteColor  = str2rgb('yellow');
    
    for i = 1:numel(MinPhaseNames)
        pname = MinPhaseNames{i};
        if strcmpi(pname,'Forsterite'), ebsd_raw(pname).color = ForsteriteColor;
        elseif strcmpi(pname,'Diopside'), ebsd_raw(pname).color = DiopsideColor;
        elseif strcmpi(pname,'Anorthite'), ebsd_raw(pname).color = AnorthiteColor;
        end
    end
    clear AnorthiteColor DiopsideColor ForsteriteColor
    
    %% --- 3.3. Checkpoint File Paths ---
    madFile   = fullfile(checkpointDir, sprintf('%s_ebsd_mad.mat', sampleName));
    cropFile  = fullfile(checkpointDir, sprintf('%s_crop.mat', sampleName));
    phaseFile = fullfile(checkpointDir, sprintf('%s_ebsd_phase.mat', sampleName));
    oriFile   = fullfile(checkpointDir, sprintf('%s_ebsd_ori.mat', sampleName));
    fillFile  = fullfile(checkpointDir, sprintf('%s_ebsd_fill.mat', sampleName));
    proFile   = fullfile(checkpointDir, sprintf('%s_ebsd_pro.mat', sampleName));
    paramFile = fullfile(checkpointDir, sprintf('%s_params.mat', sampleName));
    
    % Start with the raw EBSD data
    ebsd = ebsd_raw;
    
    %% --- 3.4. Initial Plots ---
    if runStart
        showPhaseStats(ebsd_raw, phases, 'Phase distribution before cleaning');
        plotPhaseMap(ebsd_raw, sampleName, exportPath, 'raw', params.exportRes);
        plotIPFMapPhases(ebsd_raw, sampleName, exportPath, 'raw', params.exportRes);
    end
    
    %% --- 3.5. MAD Filter ---
    fprintf('--- MAD Filter Processing ---\n');
    if runMAD
        [ebsd_mad, ~] = doMADFilter(ebsd, sampleName, exportPath);
        save(madFile,'ebsd_mad');
        fprintf('✔ Checkpoint MAD saved\n');
        ebsd = ebsd_mad;
    elseif exist(madFile,'file')
        load(madFile,'ebsd_mad'); ebsd = ebsd_mad;
        fprintf('✔ MAD checkpoint loaded.\n');
    else
        disp('MAD filter skipped and no previous file exists.');
    end
    
    %% --- 3.6. Create Sample Mask (Cropping) ---
    if runCrop
        [Nrow, Ncol] = ebsd.size;
        phases = ebsd.mineralList;
        notIndexedId = find(strcmpi(phases, 'notIndexed')); 
        indexedMask = reshape(ebsd.phaseId ~= notIndexedId, Nrow, Ncol); 
        sampleMask = false(Nrow, Ncol);
        
        % Row-wise pass
        for r = 1:Nrow
            idx = find(indexedMask(r, :));
            if isempty(idx), continue; end
            sampleMask(r, idx(1):idx(end)) = true;
        end
        % Column-wise pass
        for c = 1:Ncol
            idx = find(indexedMask(:, c));
            if isempty(idx), continue; end
            sampleMask(idx(1):idx(end), c) = true;
        end
        save(cropFile,'sampleMask');
        fprintf('✔ Sample mask created and saved.\n');
        phases = ebsd.mineralList;
        showPhaseStats(ebsd(sampleMask), phases, 'Phase distribution after cropping');
    elseif exist(cropFile,'file')
        load(cropFile,'sampleMask');
        fprintf('✔ Sample mask loaded from checkpoint.\n');
    else
        disp('Cropping skipped and no previous file exists.');
        [Nrow, Ncol] = ebsd.size;
        sampleMask = true(Nrow, Ncol);  % fallback: include all pixels
    end
    
    %% --- 3.7. Data Quality Assessment (Strict vs. Relaxed Protocol) ---
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    % Use only the masked region for assessing data quality
    validPhaseIds = ebsd.phaseId(sampleMask);
    fracNotIndexed = sum(validPhaseIds == notIndexedId) / numel(validPhaseIds);
    
    fprintf('\n--- Data Quality Assessment ---\n');
    if fracNotIndexed < 0.75
        runStrict = true;
        fprintf('✔ Data sufficiently indexed (%.2f%% notIndexed). **STRICT mode** selected.\n', fracNotIndexed*100);
    else
        runStrict = false;
        fprintf('⚠ Sparse data (%.2f%% notIndexed). **RELAXED mode** selected.\n', fracNotIndexed*100);
    end
    
    %% --- 3.8. Wild Spike Removal (WSR) ---
    fprintf('\n--- WSR Processing ---\n');
    
    % Initialise variables required for checkpoints/WSR
    ebsd_phase = []; ebsd_ori = []; phaseMapClean = []; oriQuatClean = []; protectedMask = [];
    
    if runStrict
        % Phase WSR (STRICT)
        if runPhaseWSR
            fprintf('--- Running Phase WSR (Strict) ---\n');
            [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR(ebsd, sampleName, exportPath, sampleMask);
            save(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            fprintf('✔ Phase WSR saved\n'); ebsd = ebsd_phase;
        elseif exist(phaseFile,'file')
            load(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask'); ebsd = ebsd_phase;
            fprintf('✔ Loaded existing Phase WSR checkpoint.\n');
        end
        
        % Orientation WSR
        if runOriWSR
            if isempty(ebsd_phase) && ~exist(phaseFile,'file')
                fprintf('⚠ Orientation WSR cannot run without Phase WSR. Skipping.\n');
            else
                fprintf('--- Running Orientation WSR ---\n');
                [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath, sampleMask);
                save(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean'); ebsd = ebsd_ori;
                fprintf('✔ Orientation WSR saved\n');
            end
        elseif exist(oriFile,'file')
            load(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean'); ebsd = ebsd_ori;
            fprintf('✔ Loaded existing Orientation WSR checkpoint.\n');
        end
    else % RELAXED Protocol (Aggressive Phase WSR)
        % Phase WSR (RELAXED/Aggressive)
        if runPhaseWSR
            fprintf('--- Running Phase WSR (Aggressive) ---\n');
            [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR_aggressive(ebsd, sampleName, exportPath, sampleMask);
            save(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            fprintf('✔ Phase WSR saved\n'); 
            ebsd = ebsd_phase;
        elseif exist(phaseFile,'file')
            load(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask'); ebsd = ebsd_phase;
            fprintf('✔ Loaded existing Phase WSR checkpoint.\n');
        end
        
        % Orientation WSR
        if runOriWSR
            if isempty(ebsd_phase) && ~exist(phaseFile,'file')
                fprintf('⚠ Orientation WSR cannot run without Phase WSR. Skipping.\n');
            else
                fprintf('--- Running Orientation WSR ---\n');
                [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath, sampleMask);
                save(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean'); ebsd = ebsd_ori;
                fprintf('✔ Orientation WSR saved\n');
            end
        elseif exist(oriFile,'file')
            load(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean'); ebsd = ebsd_ori;
            fprintf('✔ Loaded existing Orientation WSR checkpoint.\n');
        end
    end
    
    %% --- 3.9. Hole Filling (Unprotected Pixels) ---
    fprintf('\n--- Filling Unprotected Holes ---\n');
    if runHoleFill
        ebsd_fill = ebsd;
        if runStrict
            [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingBFS(ebsd_fill, oriQuatClean, phaseMapClean, params.radius_fill, protectedMask, sampleMask);
        else
            [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingMPF(ebsd_fill, oriQuatClean, phaseMapClean, params.radius_fill, protectedMask, sampleMask);
        end
        
        % The hole-filling functions modify phase/ori maps; wrap back to EBSD.
        ebsd_fill = EBSD(ebsd_fill, 'convert');
        showPhaseStats(ebsd_fill, phases, 'Phase distribution after hole fill');
        plotPhaseMap(ebsd_fill, sampleName, exportPath, 'unproFill', params.exportRes);
        plotIPFMapPhases(ebsd_fill, sampleName, exportPath, 'unproFill', params.exportRes);
        save(fillFile,'ebsd_fill','phaseMapClean','oriQuatClean','protectedMask'); ebsd = ebsd_fill;
        fprintf('✔ Checkpoint Hole Fill saved\n');
    elseif exist(fillFile,'file')
        load(fillFile,'ebsd_fill','phaseMapClean','oriQuatClean','protectedMask'); ebsd = ebsd_fill;
        fprintf('✔ Loaded existing HoleFilling checkpoint.\n');
    else
        disp('Hole Filling skipped and no previous file exists.');
    end
    
    %% --- 3.10. Protected Pixel Filling ---
    % Fill pixels previously set to notIndexed by the MAD/WSR filters
    if runProFill
        fprintf('\n--- Filling real grain protected pixels ---\n');
        % Initialisation is required if checkpoints were loaded and ebsd is not the current WSR output
        if isempty(phaseMapClean) || isempty(oriQuatClean)
            [phaseMapClean, oriQuatClean, protectedMask] = residuals(ebsd);
        end
        
        [ebsd_pro, phaseMapClean, oriQuatClean] = doProtectedFilling(ebsd, oriQuatClean, phaseMapClean, protectedMask, sampleMask);
        ebsd = ebsd_pro;
        showPhaseStats(ebsd, phases, 'Phase distribution after Protected Fill');
        plotPhaseMap(ebsd, sampleName, exportPath, 'proFill', params.exportRes);
        plotIPFMapPhases(ebsd, sampleName, exportPath, 'proFill', params.exportRes);
        save(proFile,'ebsd_pro','phaseMapClean','oriQuatClean','protectedMask'); 
        fprintf('✔ Checkpoint Protected Fill saved to proFile \n');
        
    elseif exist(proFile,'file')
        load(proFile,'ebsd_pro','phaseMapClean','oriQuatClean','protectedMask'); 
        ebsd = ebsd_pro;
        fprintf('✔ Loaded existing Protected Pixels Filling checkpoint.\n');
    else
        disp('Protected Pixels Filling skipped and no previous file exists.');
    end
   
    %% --- 3.11. Export Final EBSD & Parameters ---
    if runSaveFile
        outFile = fullfile(dataDir,[sampleName '_clean.ctf']);
        ebsd.export(outFile);
        fprintf('✔ Saved cleaned EBSD: %s\n', outFile);
        save(paramFile,'params'); 
        fprintf('✔ Parameters saved in a mat file\n');
    end
    
    diary off % Stop logging for the current sample
end
% end of loop over samples

%% EBSD Data Cleaning Helper Functions
% Contents:
%   1. MAD Filter
%   2. Phase WSR (Strict & Aggressive)
%   3. Orientation WSR
%   4. Mean Orientation Calculations
%   5. Hole Filling (BFS & MPF)
%   6. Protected Filling
%   7. General Utilities (Stats, Plots, Residuals)
%
% NOTE: These functions rely on the existence of 'global params' defined 
% in the calling script.
%

%% =========================================================================
%   1. MAD Filter
% =========================================================================

function [ebsd_mad, badPixels] = doMADFilter(ebsd, sampleName, exportPath)
% DOMADFILTER Applies a Mean Angular Deviation (MAD) filter.
    global params;
    fprintf('\n--- Applying MAD Filter (Threshold = %.2f rad) ---\n', params.madThreshold);
    % --- Identify bad pixels ---
    badPixels = ebsd.mad > params.madThreshold;
    numBad = sum(badPixels,'all');
    fprintf('Found %d pixels exceeding MAD threshold.\n', numBad);
    % --- Set bad pixels to notIndexed ---
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    ebsd(badPixels).phaseId = notIndexedId;
    ebsd_mad = ebsd;
    fprintf('✔ MAD filter applied: %d pixels set to notIndexed.\n', numBad);
    % --- Show stats and plots ---
    showPhaseStats(ebsd_mad, phases, 'Phase distribution after MAD filter');
    plotPhaseMap(ebsd_mad, sampleName, exportPath, 'MADfilter', params.exportRes);
    plotIPFMapPhases(ebsd_mad, sampleName, exportPath, 'MADfilter', params.exportRes);
end

%% =========================================================================
%   2. Phase Wild Spike Removal (WSR)
% =========================================================================

function [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR(ebsd, sampleName, exportPath, sampleMask)
% DOPHASESWSR Performs Phase WSR (Strict Protocol).
    global params;
    fprintf('Starting Phase WSR (Radius = %d) \n', params.radius_phase);
    phases = ebsd.mineralList;
    Nrow = ebsd.size(1); Ncol = ebsd.size(2);
    phaseMapoG = reshape(double(ebsd.phaseId), Nrow, Ncol);
    phaseMapClean = phaseMapoG;
    % Identify notIndexed pixels
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    % Build quaternion grid
    oriQuatClean = zeros(Nrow,Ncol,4);
    N_total = Nrow*Ncol;
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); 
        mask = (phaseMapoG == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        q = quaternion(ebsd(mask).orientations);
        idx = find(mask);
        oriQuatClean(idx)             = [q.a];  
        oriQuatClean(idx + N_total)   = [q.b];
        oriQuatClean(idx + 2*N_total) = [q.c]; 
        oriQuatClean(idx + 3*N_total) = [q.d];
    end
    % Precompute kernel
    kernel_phase = double(fspecial('disk', params.radius_phase))>0;
    kernel_phase(params.radius_phase+1, params.radius_phase+1) = 0;
    protectedMask = ebsd.mad > params.madThreshold; % Keep MAD-filtered pixels protected
    % Loop over phase
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        fprintf('Processing phase: %s\n', pname);
        mask = (phaseMapoG == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        [rows, cols] = find(mask);
        paddedPhase = padarray(phaseMapClean,[params.radius_phase params.radius_phase],0,'both');
        paddedOri   = padarray(oriQuatClean,[params.radius_phase params.radius_phase],0,'both');
        % loop over pixels
        for k = 1:numel(rows)
            i = rows(k); j = cols(k);
            iP = i + params.radius_phase; jP = j + params.radius_phase;
            win = paddedPhase(iP-params.radius_phase:iP+params.radius_phase, ...
                              jP-params.radius_phase:jP+params.radius_phase);
            neigh = win(kernel_phase);
            total_neigh = numel(neigh);
            indexedMask = (neigh > 0 & neigh ~= notIndexedId);
            validNeigh = neigh(indexedMask);
            Ni = numel(validNeigh)/total_neigh;
            % CASE 1: few valid neighbours → set non-indexed
            if Ni <= 0.25
                phaseMapClean(i,j) = notIndexedId;
                oriQuatClean(i,j,:) = 0;
                protectedMask(i,j) = true;
                continue;
            end
            % CASE 2: 25% < Ni < 50%
            if Ni > 0.25 && Ni < 0.5
                uniquePhases = unique(validNeigh);
                if isscalar(uniquePhases)
                    maj = uniquePhases;
                    if phaseMapClean(i,j) ~= maj
                        phaseMapClean(i,j) = maj;
                    end
                    % Update orientation
                    oriWin = paddedOri(iP-params.radius_phase:iP+params.radius_phase, ...
                                       jP-params.radius_phase:jP+params.radius_phase,:);
                    oriList = reshape(oriWin, [], 4);
                    validOriMask = (win(kernel_phase) == maj);
                    neighbourQuatsList = oriList(validOriMask,:);
                    Nneighbours = size(neighbourQuatsList,1);
                    if Nneighbours >= 2
                        currentQ_vec = squeeze(oriQuatClean(i,j,:))';
                        [meanOri, ~] = calc_mean_ori_wsr_normal(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                        oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                    end
                end
                continue;
            end
            % CASE 3: spike fix based on dominant phase
            fracThresh = (Ni >= 0.75) * (2/3) + (Ni < 0.75) * (3/4);
            maj = mode(validNeigh);
            fracMaj = sum(validNeigh == maj)/numel(validNeigh);
            if maj ~= phaseMapClean(i,j) && fracMaj >= fracThresh
                phaseMapClean(i,j) = maj;
                oriWin = paddedOri(iP-params.radius_phase:iP+params.radius_phase, ...
                                   jP-params.radius_phase:jP+params.radius_phase,:);
                oriList = reshape(oriWin, [], 4);
                validOriMask = (win(kernel_phase) == maj);
                neighbourQuatsList = oriList(validOriMask,:);
                Nneighbours = size(neighbourQuatsList,1);
                if Nneighbours >= 2
                    currentQ_vec = squeeze(oriQuatClean(i,j,:))';
                    [meanOri, ~] = calc_mean_ori_wsr_normal(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                    oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                end
            end
        end % end loop over pixels
    end % end loop over phases
    % Update EBSD object
    ebsd_phase = ebsd;
    ebsd_phase.phaseId(:) = phaseMapClean(:);
    qFull_phase = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        mask = (phaseMapClean==pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        ebsd_phase(mask).orientations = orientation(qFull_phase(mask), ebsd(pname).CS);
    end
    % Show stats and plots
    showPhaseStats(ebsd_phase, phases, 'Phase distribution after Phase WSR');
    plotPhaseMap(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    plotIPFMapPhases(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);   
end

function [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR_aggressive(ebsd, sampleName, exportPath, sampleMask)
% DOPHASESWSR_AGGRESSIVE Performs Phase WSR (Aggressive/Relaxed Protocol).
    global params;
    minDomFrac = params.min_dom_frac;
      
    fprintf('\n--- Starting Phase WSR (Radius = %d, Aggressive) ---\n', params.radius_phase);
    phases = ebsd.mineralList;
    Nrow = ebsd.size(1); Ncol = ebsd.size(2);
    phaseMapoG = reshape(double(ebsd.phaseId), Nrow, Ncol);
    phaseMapClean = phaseMapoG;
    
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    
    % Build quaternion grid 
    oriQuatClean = zeros(Nrow,Ncol,4);
    N_total = Nrow*Ncol;
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); 
        mask = (phaseMapoG == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        q = quaternion(ebsd(mask).orientations);
        idx = find(mask);
        oriQuatClean(idx)             = [q.a];  
        oriQuatClean(idx + N_total)   = [q.b];
        oriQuatClean(idx + 2*N_total) = [q.c]; 
        oriQuatClean(idx + 3*N_total) = [q.d];
    end
    
    % Precompute kernel
    kernel_phase = double(fspecial('disk', params.radius_phase))>0;
    kernel_phase(params.radius_phase+1, params.radius_phase+1) = 0;
    % protect hihg_MAD pixels
    protectedMask = ebsd.mad > params.madThreshold; 
    
    % Loop over phase
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        fprintf('Processing phase: %s\n', pname);
        mask = (phaseMapoG == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
       
        [rows, cols] = find(mask);
        paddedPhase = padarray(phaseMapClean,[params.radius_phase params.radius_phase],0,'both');
        paddedOri   = padarray(oriQuatClean,[params.radius_phase params.radius_phase],0,'both');
        
        countRemoved = 0;
        countPhaseFlip = 0;
        
        % loop over pixels
        for k = 1:numel(rows)
            i = rows(k); j = cols(k);
            iP = i + params.radius_phase; jP = j + params.radius_phase;
            win = paddedPhase(iP-params.radius_phase:iP+params.radius_phase, ...
                              jP-params.radius_phase:jP+params.radius_phase);
            neigh = win(kernel_phase);
            total_neigh = numel(neigh);
            indexedMask = (neigh > 0 & neigh ~= notIndexedId);
            validNeigh = neigh(indexedMask);
            Ni = numel(validNeigh)/total_neigh;
            
            % --- CASE 1: MINIMAL REMOVAL CHECK ---
            if Ni <= 0.10
                phaseMapClean(i,j) = notIndexedId;
                oriQuatClean(i,j,:) = 0;
                protectedMask(i,j) = true; 
                countRemoved = countRemoved + 1;
                continue; % Must skip if data is invalid
            end
            
            % --- CASE 2: AGGRESSIVE PHASE SPIKE FIX ---
            [uniquePh, ~, ic] = unique(validNeigh);
            counts = accumarray(ic,1);
            [maxCount, idxMax] = max(counts);
            domPhase = uniquePh(idxMax);
            fracMaj = maxCount / numel(validNeigh);
            
            if phaseMapClean(i,j) ~= domPhase && fracMaj >= minDomFrac
                phaseMapClean(i,j) = domPhase;
                countPhaseFlip = countPhaseFlip + 1;
                
                % Update orientation based on new phase
                oriWin = paddedOri(iP-params.radius_phase:iP+params.radius_phase, ...
                                   jP-params.radius_phase:jP+params.radius_phase,:);
                oriList = reshape(oriWin, [], 4);
                validOriMask = (win(kernel_phase) == domPhase);
                neighbourQuatsList = oriList(validOriMask,:);
                Nneighbours = size(neighbourQuatsList,1);
                
                if Nneighbours >= 2
                    currentQ_vec = squeeze(oriQuatClean(i,j,:))';
                    % Must call the aggressive Lead-Check version for smoothing
                    [meanOri, ~] = calc_mean_ori_wsr_aggressive(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                    oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                end
            end
        end % End loop over pixels
        
        % Print summary statement
        fprintf('  ✔ Phase %s: Removed %d points. Flipped %d points.\n', pname, countRemoved, countPhaseFlip);
        
    end % End loop over phases
    
    % Update EBSD object
    ebsd_phase = ebsd;
    ebsd_phase.phaseId(:) = phaseMapClean(:);
    qFull_phase = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        mask = (phaseMapClean==pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        ebsd_phase(mask).orientations = orientation(qFull_phase(mask), ebsd(pname).CS);
    end
    
    % Show stats and plots
    plotPhaseMap(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    plotIPFMapPhases(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
end

%% =========================================================================
%   3. Orientation Wild Spike Removal (WSR)
% =========================================================================

function [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath, sampleMask)  
% DOORITENTATIONWSR Performs Orientation WSR, including twin checks for Anorthite.
    global params;
    radius_ori  = params.radius_ori;
    misTol_ori  = params.misTol_ori;
    minFrac_ori = params.minFrac_ori;
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol] = size(phaseMapClean);
    kernel_ori = double(fspecial('disk', radius_ori)) > 0;
    kernel_ori(radius_ori+1, radius_ori+1) = 0;
    % loop over phase
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid) & sampleMask;
        [rowsAll, colsAll] = find(mask);
        numPhasePixels = numel(rowsAll);
        % Pad orientation and phase maps for neighbourhood access
        paddedOri   = padarray(oriQuatClean, [radius_ori radius_ori], 0, 'both');
        paddedPhase = padarray(phaseMapClean, [radius_ori radius_ori], 0, 'both');
        % --- Wild spike pre-filter ---
        wildSpikes = false(Nrow, Ncol);
        for k = 1:numPhasePixels
            i = rowsAll(k); j = colsAll(k);
            iP = i + radius_ori; jP = j + radius_ori;
            oriWin   = paddedOri(iP-radius_ori:iP+radius_ori, jP-radius_ori:jP+radius_ori, :);
            phaseWin = paddedPhase(iP-radius_ori:iP+radius_ori, jP-radius_ori:jP+radius_ori);
            neighbourMask = (phaseWin == pid) & kernel_ori;
            oriList = reshape(oriWin, [], 4);
            neighbourQuatsList = oriList(neighbourMask(:), :);
            if isempty(neighbourQuatsList)
                continue;
            end  
            currentQ_vec = squeeze(oriQuatClean(i,j,:))';
            dots = abs(neighbourQuatsList * currentQ_vec'); dots(dots>1)=1;
            misAngles = 2*acos(dots);
            fracSimilar = sum(misAngles < misTol_ori) / numel(misAngles);
            if fracSimilar < minFrac_ori
                wildSpikes(i,j) = true;
            end
        end
        % --- Process only wild spike pixels ---
        [rowsSpike, colsSpike] = find(wildSpikes);
        numWildSpikes = numel(rowsSpike);
        fprintf('Phase %s: total pixels = %d, potential wild spikes = %d\n', pname, numPhasePixels, numWildSpikes);
        % initialise counters
        numRemoved    = 0;
        numReoriented = 0;
        numSkippedTwin = 0;
        for k = 1:numWildSpikes
            i = rowsSpike(k); j = colsSpike(k);
            iP = i + radius_ori; jP = j + radius_ori;
            oriWin   = paddedOri(iP-radius_ori:iP+radius_ori, jP-radius_ori:jP+radius_ori, :);
            phaseWin = paddedPhase(iP-radius_ori:iP+radius_ori, jP-radius_ori:jP+radius_ori);
            neighbourMask = (phaseWin == pid) & kernel_ori;
            oriList = reshape(oriWin, [], 4);
            neighbourQuatsList = oriList(neighbourMask(:), :);
            Nneighbours = size(neighbourQuatsList,1);
            currentQ_vec = squeeze(oriQuatClean(i,j,:))';
            currentQ = quaternion(currentQ_vec(1), currentQ_vec(2), currentQ_vec(3), currentQ_vec(4));
            % --- 0 neighbours ---
            if Nneighbours == 0
                phaseMapClean(i,j) = notIndexedId;
                oriQuatClean(i,j,:) = 0;
                paddedOri(iP,jP,:)  = 0;
                numRemoved = numRemoved + 1;
                continue;
            end
            % --- 1 neighbour ---
            if Nneighbours == 1, continue; end
            % --- 2 neighbours ---
            if Nneighbours == 2
                dots = abs(neighbourQuatsList * currentQ_vec'); dots(dots>1)=1;
                mori = 2*acos(dots);
                if all(mori < misTol_ori)
                    qMean = mean(quaternion(neighbourQuatsList(:,1), neighbourQuatsList(:,2), ...
                                            neighbourQuatsList(:,3), neighbourQuatsList(:,4)), 'meanOrientation');
                    oriQuatClean(i,j,:) = [qMean.a qMean.b qMean.c qMean.d];
                    paddedOri(iP,jP,:)  = [qMean.a qMean.b qMean.c qMean.d];
                    numReoriented = numReoriented + 1;
                end
                continue;
            end
            % --- N >= 3 neighbours ---
            dots = abs(neighbourQuatsList * currentQ_vec'); dots(dots>1)=1;
            misAngles = 2*acos(dots);
            if all(misAngles < misTol_ori)
                qMean = mean(quaternion(neighbourQuatsList.'), 'meanOrientation');
            else
                qMean = calc_mean_ori_wsr_normal(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
            end
            % --- Twin check for Anorthite ---
            if strcmpi(pname,'Anorthite')
                % Define twin laws (you can move these outside the loop for efficiency)
                cs = ebsd('Anorthite').CS;
                twinLaws = { ...
                    {'Albite',    orientation.byAxisAngle(vector3d(0,1,0), 180*degree, cs), 5*degree}, ...
                    {'Pericline', orientation.byAxisAngle(vector3d(1,0,0), 180*degree, cs), 5*degree}, ...
                    {'Carlsbad',  orientation.byAxisAngle(vector3d(0,0,1), 180*degree, cs), 5*degree}, ...
                    {'Manebach',  orientation(reflection(Miller(0,0,1,cs))), 5*degree}, ...
                    {'Baveno',    orientation(reflection(Miller(0,2,1,cs))), 5*degree} ...
                };
                % Compute misorientation between current and mean
                misOri = qMean * inv(currentQ);
                % Check if misorientation is close to any twin law
                isTwin = false;
                for t = 1:numel(twinLaws)
                    law = twinLaws{t};
                    if angle(misOri, law{2}) < law{3}, isTwin = true; break; end
                end
                if isTwin, numSkippedTwin = numSkippedTwin + 1; continue; end % Skip if it’s a twin
            end
            % --- Update orientation ---
            oriQuatClean(i,j,:) = [qMean.a qMean.b qMean.c qMean.d];
            paddedOri(iP,jP,:)  = [qMean.a qMean.b qMean.c qMean.d];
            numReoriented = numReoriented + 1;
        end % end loop over wild spikes
        % --- After loop, print summary ---
        fprintf('Phase %s Ori WSR summary: Removed = %d, Reoriented = %d, Skipped twins = %d\n', ...
        pname, numRemoved, numReoriented, numSkippedTwin);
    end % end loop over phases
    % --- Push back to EBSD ---
    ebsd_ori = ebsd;
    ebsd_ori.phaseId(:) = phaseMapClean(:);
    qFull_wsr = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        mask = (phaseMapClean==pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        ebsd_ori(mask).orientations = orientation(qFull_wsr(mask), ebsd_ori(pname).CS);
    end
    % Show stats and plots
    showPhaseStats(ebsd_ori, phases, 'Phase distribution after Ori WSR');
    plotPhaseMap(ebsd_ori, sampleName, exportPath, 'OriWSR', params.exportRes);
    plotIPFMapPhases(ebsd_ori, sampleName, exportPath, 'OriWSR', params.exportRes);
end

%% =========================================================================
%   4. Mean Orientation Calculations
% =========================================================================

function [meanOri, clusterSizes] = calc_mean_ori_wsr_normal(qList, misTol_ori, Nneighbours, currentQ_vec)
% CALC_MEAN_ORI_WSR_NORMAL Calculates mean orientation using hierarchical clustering and threshold fraction (Strict).
    global params
    thresholdFrac = params.thresholdFrac;
    % Convert current orientation vector to quaternion
    currentQ = quaternion(currentQ_vec(1), currentQ_vec(2), ...
                          currentQ_vec(3), currentQ_vec(4));
                          
    % --- Step 1: pairwise angular distances ---
    D = 2 * acos(min(abs(qList*qList.'),1)); % radians
    D(1:size(D,1)+1:end) = 0;                % diagonal 0 manually
    Dc = squareform(D);                      % condensed distance matrix
    
    % --- Step 2: hierarchical clustering ---
    Z = linkage(Dc,'single');
    
    % --- Step 3: cluster neighbours by misorientation tolerance ---
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');
    
    % --- Step 4: select dominant cluster ---
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;
    
    % Step 5: compute representative orientation based on dominant cluster fraction
    domClusterFrac = maxCount / Nneighbours;
    
    % --- CASE 1: dominant cluster is sufficiently strong ---
    if domClusterFrac >= thresholdFrac
        qCluster = quaternion(qList(members,:).'); % Transpose Mx4 for MTEX
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    end
    
    % --- CASE 2: dominant cluster is weak → pick cluster mean closest to current pixel ---
    uniqueClusters = unique(idx);
    minMisorientation = inf;
    closestClusterMean = quaternion(0,0,0,1); % placeholder
    for c = 1:numel(uniqueClusters)
        clusterMembers = (idx == uniqueClusters(c));
        qCluster = quaternion(qList(clusterMembers,:).'); % 4xM
        clusterMean = mean(qCluster, 'meanOrientation');
        % misorientation with current pixel
        mis = angle(clusterMean * currentQ);
        if mis < minMisorientation
            minMisorientation = mis;
            closestClusterMean = clusterMean;
        end
    end
    meanOri = closestClusterMean;
end

function [meanOri, clusterSizes] = calc_mean_ori_wsr_aggressive(qList, misTol_ori, Nneighbours, currentQ_vec)
% CALC_MEAN_ORI_WSR_AGGRESSIVE Calculates mean orientation using hierarchical clustering and Lead check (Aggressive).
    global params
    minLead = params.minLead;
    scaleLead = params.scaleLead;                         
    % Convert current orientation vector to quaternion
    currentQ = quaternion(currentQ_vec(1), currentQ_vec(2), ...
                          currentQ_vec(3), currentQ_vec(4));        
    % --- Step 1 & 2: Pairwise distances and hierarchical clustering (Unchanged) ---
    D = 2 * acos(min(abs(qList*qList.'),1)); % radians
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);
    Z = linkage(Dc,'single');
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');
    % --- Step 3: Dominant cluster check (Lead) ---
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;
    % Calculate Lead
    otherCounts = counts; otherCounts(domCluster) = 0;
    leadDiff = maxCount - max(otherCounts, 0); % Absolute difference  
    % Calculate Required Lead (Relaxed/Scaled requirement, same as hole filling)
    reqLead = max(minLead, ceil(scaleLead * Nneighbours));
    % --- Step 4: Success check using Lead  ---
    if leadDiff >= reqLead
        % CASE 1: Dominant cluster is sufficiently strong by Lead
        qCluster = quaternion(qList(members,:).'); % Transpose Mx4 for MTEX
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    else
        % CASE 2: Dominant cluster is weak (Lead Check failed). 
        uniqueClusters = unique(idx);
        minMisorientation = inf;
        closestClusterMean = quaternion(0,0,0,1); % Placeholder     
        for c = 1:numel(uniqueClusters)
            clusterMembers = (idx == uniqueClusters(c));
            qCluster = quaternion(qList(clusterMembers,:).');
            clusterMean = mean(qCluster, 'meanOrientation');            
            % misorientation with current pixel
            mis = angle(clusterMean * currentQ);
            if mis < minMisorientation
                minMisorientation = mis;
                closestClusterMean = clusterMean;
            end
        end
        meanOri = closestClusterMean;
        return;
    end
end

%% =========================================================================
%   5. Hole Filling (Unprotected)
% =========================================================================

function [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingBFS(ebsd, oriQuatClean, phaseMapClean, radii, protectedMask, sampleMask)
% DOHOLEFILLINGBFS Performs Breadth-First Search (BFS) based hole filling (Strict Protocol).
    global params;
    misTol_ori = params.misTol_ori;
    phaseFrac_all = params.phaseFrac;   % map(radius) -> [Ni_thresh, fracDom_thresh]
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol] = size(phaseMapClean);
    
    % Loop over radii 
    for radius = radii
        fprintf('--- Hole Filling: radius = %d ---\n', radius);
        phaseFrac   = phaseFrac_all(radius);  % [Ni_threshold, fracDom_threshold]
        
        % Precompute full disk kernel and outerKernel (exclude centre)
        N = 2*radius + 1;
        kernelFull = double(fspecial('disk', radius)) > 0;
        innerKernel = false(N); innerKernel(radius+1, radius+1) = true;
        outerKernel = kernelFull & ~innerKernel;
        numNeighbours = sum(outerKernel(:));   % number of outer neighbours in full kernel
        
        % Per-radius bookkeeping
        skipMask = false(Nrow, Ncol);  % true => this pixel failed fill for this radius
        visited  = false(Nrow, Ncol);  % true => this pixel has been discovered as part of a cluster already
        holeMask = (phaseMapClean == notIndexedId) & ~protectedMask & sampleMask;
        fillCountRadius = 0;   % track successful fills (orientation+phase)
        clusterId = 0;
        
        % While there exists an undiscovered hole, form a cluster and fill it
        while true
            % find first undiscovered hole
            [r0, c0] = find(holeMask & ~visited, 1, 'first');
            if isempty(r0)
                fprintf('⚠ No more undiscovered holes in radius = %d \n', radius);
                 fprintf('✔ Total filled in this radius: %d\n', fillCountRadius);
                break; % no more undiscovered holes at this radius
            end
            clusterId = clusterId + 1;
           
            % ---------- DISCOVERY BFS (gather connected component) ----------
            % Use 8-connectivity for cluster definition
            queue = [r0, c0];
            visited(r0,c0) = true;
            head = 1;
            while head <= size(queue,1)
                i = queue(head,1); j = queue(head,2);
                head = head + 1;
                % explore 8 neighbours
                for di = -1:1
                    for dj = -1:1
                        if di==0 && dj==0, continue; end
                        ni = i + di; nj = j + dj;
                        if ni<1 || nj<1 || ni>Nrow || nj>Ncol, continue; end
                        % neighbour is a hole (notIndexed) and not protected
                        if holeMask(ni,nj) && ~visited(ni,nj)
                            visited(ni,nj) = true;
                            queue(end+1, :) = [ni, nj];
                        end
                    end
                end
            end
            % ---------- FILL PHASE ----------
            totalCandidates = size(queue,1);
            filledCount = 0;
            % Iterate cluster pixels
            for q = 1:totalCandidates
                i = queue(q,1); j = queue(q,2);
                % If pixel already filled by earlier operation in this radius, skip
                if phaseMapClean(i,j) ~= notIndexedId
                    continue;
                end
                % If previously marked as skip for this radius, skip attempt
                if skipMask(i,j)
                    continue;
                end
                % Build local window bounds (clamp at image edges)
                rmin = max(i - radius, 1); rmax = min(i + radius, Nrow);
                cmin = max(j - radius, 1); cmax = min(j + radius, Ncol);
                winPhase = phaseMapClean(rmin:rmax, cmin:cmax);
                winOri   = oriQuatClean(rmin:rmax, cmin:cmax, :);
                % Valid neighbour mask inside window: indexed and not 'notIndexed'
                validMask = (winPhase > 0 & winPhase ~= notIndexedId) & outerKernel(1:(rmax-rmin+1),1:(cmax-cmin+1));
                neighPhases = winPhase(validMask);
                
                % Condition 1: check fraction of valid neighbours (Ni threshold)
                Ni = numel(neighPhases) / max(1, numNeighbours); 
                if Ni < phaseFrac(1)
                    skipMask(i,j) = true; % cannot fill with this radius
                    continue;
                end
                % Condition 2: dominant phase fraction among valid neighbours
                [uniquePh, ~, ic] = unique(neighPhases);
                counts = accumarray(ic, 1);
                [maxCount, idxMax] = max(counts);
                domPhase = uniquePh(idxMax);
                fracDom = maxCount / max(1, numel(neighPhases));
                if fracDom < phaseFrac(2)
                    skipMask(i,j) = true;
                    continue;
                end
                
                % Compute mean orientation using q's of valid dominated-phase neighbours
                qList = reshape(winOri, [], 4); % rows correspond to same ordering as winPhase(:)
                qMean = calc_mean_ori_hole_strict(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, validMask);
                 if isempty(qMean)
                    skipMask(i,j) = true;
                    continue;
                end
      
                % ===== SUCCESSFUL FILL =====
                phaseMapClean(i,j) = domPhase;
                oriQuatClean(i,j,:) = [qMean.a qMean.b qMean.c qMean.d];
                filledCount     = filledCount + 1;
                fillCountRadius = fillCountRadius   + 1;
            end
            % Print cluster result if filled with anything
            if filledCount > 0 && totalCandidates > 10
                fprintf('Cluster %d: filled %d/%d\n', clusterId, filledCount, totalCandidates);
            end
            % Update holeMask for next cluster (we exclude protected pixels and those already filled)
            holeMask = (phaseMapClean == notIndexedId) & ~protectedMask & sampleMask;
        end % end while true (BFS cluster loop)
    end % end loop over radii
    
    % Rebuild EBSD object with final maps
    ebsd_fill = ebsd;
    ebsd_fill.phaseId(:) = phaseMapClean(:);
    qFull_fill = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        ebsd_fill(mask).orientations = orientation(qFull_fill(mask), ebsd(pname).CS);
    end
end

function [meanOri, clusterSizes] = calc_mean_ori_hole_strict(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, maskValid)
% CALC_MEAN_ORI_HOLE_STRICT Calculates mean orientation for strict hole filling (Fraction Threshold + Ring Fallback).
    global params
    thresholdFrac = params.thresholdFrac;  % Fraction of dominant cluster to accept
    % Filter to dominant phase neighbours
    qListDom = qList(neighPhases == domPhase, :);
    Ndom = size(qListDom,1);
    
    if Ndom == 0
        meanOri = [];   % No valid neighbors → do not assign orientation
        return;
    end
     if Ndom == 1
        if isempty(qListDom)
            meanOri = []; % safe fallback
        else
            meanOri = quaternion(qListDom(1,1), qListDom(1,2), qListDom(1,3), qListDom(1,4));
        end
        clusterSizes = 1;
        return;
    end
    
    % --- Step 1: pairwise angular distances ---
    D = 2 * acos(min(abs(qListDom*qListDom.'), 1));
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);
    
    % --- Step 2: Hierarchical clustering by misoreintation ---
    Z = linkage(Dc,'single');
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');
    
    % --- Step 3: dominant cluster check (Fraction) ---
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx==domCluster);
    clusterSizes = counts;
    domFrac = maxCount / Ndom;
    if domFrac >= thresholdFrac
        qCluster = quaternion(qListDom(members,:).');
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    end
    
    % --- Step 4: ring-based fallback ---
    distMap = bwdist(innerKernel);
    distMap = round(distMap);
    distMap(~outerKernel) = 0;
    found = false;
    for ringWidth = 2:(radius-1)
        ringMask = (distMap == ringWidth);
        if ~any(ringMask(:)), continue; end
        % Indices in the window that are valid and in the ring
        validLinearIdx = find(maskValid);
        ringLinearIdx = find(ringMask(:));
        [~, commonidx, ~] = intersect(validLinearIdx, ringLinearIdx);         
        % Pick only dominant phase members
        qRing = qList(commonidx, :);  
        neighPhasesRing = neighPhases(commonidx);
        qRingDom = qRing(neighPhasesRing == domPhase, :);
        Nring = size(qRingDom,1);
        if Nring < 2, continue; end
        % Cluster dominant phase in the ring
        D_ring = 2*acos(min(abs(qRingDom*qRingDom.'),1));
        D_ring(1:size(D_ring,1)+1:end) = 0;
        Dc_ring = squareform(D_ring);
        Z_ring = linkage(Dc_ring,'single');
        idx_ring = cluster(Z_ring,'cutoff',misTol_ori,'criterion','distance');
        counts_ring = accumarray(idx_ring,1);
        [maxCountRing, domClusterRing] = max(counts_ring);
        % Check if dominant cluster fraction is sufficient
        domFracRing = maxCountRing / Nring;
        if domFracRing >= 0.5
            membersRing = (idx_ring==domClusterRing);  
            qClusterRing = quaternion(qRingDom(membersRing,:).'); 
            meanOri = mean(qClusterRing,'meanOrientation');
            found = true;
            break;
        end
    end
    % --- Step 5: fallback ---
    if ~found
        meanOri = [];
    end
end

function [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingMPF(ebsd, oriQuatClean, phaseMapClean, radii, protectedMask, sampleMask)
% DOHOLEFILLINGMPF Performs Multi-Pass Filler (MPF) based hole filling (Relaxed Protocol).
    global params
    scaleLead = params.scaleLead;
    minLead = params.minLead;
    misTol_ori  = params.misTol_ori;
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol] = size(phaseMapClean);
    ebsdPhaseIdGrid = reshape(ebsd.phaseId, Nrow, Ncol);
    initialHoleMask = (ebsdPhaseIdGrid == notIndexedId) & ~protectedMask & sampleMask;
    % Loop over radii 
    for radius = radii
        phaseFrac = params.phaseFrac(radius);
        N = 2 * radius + 1;
        kernelFull = double(fspecial('disk', radius)) > 0;
        innerKernel = false(N); innerKernel(radius+1, radius+1) = true;
        outerKernel = kernelFull & ~innerKernel;
        passes = 0;
        filledLastPass = false(Nrow, Ncol);
        % Multiple passes for filling holes - run until no holes
        fillCountRadius = 0;
        while true
            passes = passes + 1;
            fprintf('\n--- Hole Filling Pass %d (Radius=%d) ---\n', passes, radius);
            
            % --- Initialise padded maps for current pass ---
            % These must be updated every pass as the underlying data changes
            paddedPhase   = padarray(phaseMapClean, [radius radius], 0, 'both');
            paddedOri     = padarray(oriQuatClean, [radius radius], 0, 'both');
            % Initial hole set: All unindexed, unprotected pixels within the sample.
            baseHoleMask = (phaseMapClean == notIndexedId) & ~protectedMask & sampleMask;
            if passes == 1
                % Pass 1: Check all current holes
                holeMask = baseHoleMask;
            else
                % Subsequent Passes: Only check holes that are adjacent to the previous fills
                dilatedFilled = imdilate(filledLastPass, outerKernel);
                holeMask = baseHoleMask & dilatedFilled; % Only holes adjacent to last fills
            end
            
            [holeRows, holeCols] = find(holeMask);
            numHoles = numel(holeRows);
            % If no holes → done for this radius
            if numHoles == 0
                fprintf('⚠ No more undiscovered holes in radius = %d \n', radius); 
                fprintf('✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
                break; 
            end
            
            fillCountPass = 0;
            filledLastPass = false(Nrow, Ncol);  % Track pixels filled in this pass
            fprintf('Total candidate holes: %d\n', numHoles);
            % Loop over each hole
            for k = 1:numHoles
                i = holeRows(k); j = holeCols(k);
                if phaseMapClean(i,j) ~= notIndexedId
                    continue;   % has been filled earlier this pass or externally
                end
                
                iP = i + radius; jP = j + radius;
                % Extract local window
                winPhase = paddedPhase(iP-radius:iP+radius, jP-radius:jP+radius);
                winOri   = paddedOri(iP-radius:iP+radius, jP-radius:jP+radius, :);
                
                % Mask valid neighbours
                maskValid = (winPhase > 0 & winPhase ~= notIndexedId) & outerKernel;
                neighPhases = winPhase(maskValid);
                
                numNeighboursTotal = sum(outerKernel(:)); % Total neighbours in the kernel
                Ni = numel(neighPhases) / numNeighboursTotal;
                if Ni < phaseFrac(1), continue; end
               
                % Dominant phase
                [uniquePh, ~, ic] = unique(neighPhases);
                counts = accumarray(ic,1);
                [maxCount, idxMax] = max(counts);
                domPhase = uniquePh(idxMax);
                fracDom = maxCount / numel(neighPhases);
                if fracDom < phaseFrac(2), continue; end
               
                % Dominant Orientations
                oriList = reshape(winOri, [], 4);
                qList = oriList(maskValid(:), :);
                qMean = calc_mean_ori_hole_relaxed(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, maskValid, scaleLead, minLead);
                % Only assign phase & orientation if qMean is valid
                if ~isempty(qMean)
                    phaseMapClean(i,j) = domPhase;
                    oriQuatClean(i,j,:) = [qMean.a qMean.b qMean.c qMean.d];
                    
                    % Update padded maps immediately to allow influence in later checks within the pass
                    paddedPhase(iP,jP) = domPhase;
                    paddedOri(iP,jP,:) = [qMean.a qMean.b qMean.c qMean.d];
                    
                    filledLastPass(i,j) = true;
                    fillCountPass = fillCountPass + 1;
                    fillCountRadius = fillCountRadius + 1;
                end
            end
            fprintf('✔ Pass %d complete: filled %d holes\n', passes, fillCountPass);
            
            % Stop if no holes filled (stabilised)
            if fillCountPass == 0
                fprintf('⚠ None of the discovered holes filled.\n'); 
                fprintf('✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
                break;
            end
        end % end while (MultiPass)
    end % end loop over radii
    
    % Rebuild EBSD object
    ebsd_fill = ebsd;
    ebsd_fill.phaseId(:) = phaseMapClean(:);
    qFull_fill = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (ebsd_fill.phaseId == pid) & initialHoleMask(:);
        if ~any(mask), continue; end
        ebsd_fill(mask).orientations = orientation(qFull_fill(mask), ebsd(pname).CS);
    end
end

function [meanOri, clusterSizes] = calc_mean_ori_hole_relaxed(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, maskValid, scaleLead, minLead)
% CALC_MEAN_ORI_HOLE_RELAXED Calculates mean orientation for relaxed hole filling (Lead Check + Ring Fallback).
    clusterSizes = [];
    
    domIdx = (neighPhases == domPhase);
    qListDom = qList(domIdx, :);
    Ndom = size(qListDom,1);
    
    if Ndom == 0
        meanOri = [];
        return;
    end
    if Ndom == 1
        meanOri = quaternion(qListDom(1,1), qListDom(1,2), qListDom(1,3), qListDom(1,4));
        clusterSizes = 1;
        return;
    end
    
    % --- Step 1 & 2: Initial clustering ---
    D = 2 * acos(min(abs(qListDom*qListDom.'), 1));
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);
    Z = linkage(Dc,'single');
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');
    
     % --- Step 3: Dominant cluster check (Lead) ---
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx==domCluster);
    clusterSizes = counts;
    otherCounts = counts; otherCounts(domCluster) = 0;
    leadDiff = maxCount - max(otherCounts,0);
    
    if leadDiff >= minLead
        qCluster = quaternion(qListDom(members,1), qListDom(members,2), qListDom(members,3), qListDom(members,4));
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    end
    
    % --- Step 4: ring-based fallback ---
    distMap = bwdist(innerKernel);
    distMap = round(distMap);
    distMap(~outerKernel) = 0;
    found = false;
    for ringWidth = 2:(radius-1)
        ringMask = (distMap == ringWidth);
        if ~any(ringMask(:)), continue; end
        % Indices in the window that are valid and in the ring
        validLinearIdx = find(maskValid);
        ringLinearIdx = find(ringMask(:));
        [~, commonidx, ~] = intersect(validLinearIdx, ringLinearIdx);         
        % Pick only dominant phase members
        qRing = qList(commonidx, :);  
        neighPhasesRing = neighPhases(commonidx);
        qRingDom = qRing(neighPhasesRing == domPhase, :);
        Nring = size(qRingDom,1);
        if Nring < 2, continue; end
        % Cluster dominant phase in the ring
        D_ring = 2*acos(min(abs(qRingDom*qRingDom.'),1));
        D_ring(1:size(D_ring,1)+1:end) = 0;
        Dc_ring = squareform(D_ring);
        Z_ring = linkage(Dc_ring,'single');
        idx_ring = cluster(Z_ring,'cutoff',misTol_ori,'criterion','distance');
        counts_ring = accumarray(idx_ring,1);
        [maxCountRing, domClusterRing] = max(counts_ring);
        % Lead check in ring
        other_ring = counts_ring; other_ring(domClusterRing) = 0;
        leadDiffRing = maxCountRing - max(other_ring,0);
        reqLeadRing = max(minLead, ceil(scaleLead*Nring));
        if leadDiffRing >= reqLeadRing
            membersRing = (idx_ring == domClusterRing);
            qClusterRing = quaternion(qRingDom(membersRing,:).');
            meanOri = mean(qClusterRing,'meanOrientation');
            found = true;
            break;
        end
    end
    
    % --- Step 5: fallback ---
    if ~found
        meanOri = [];
    end
end

%% =========================================================================
%   6. Protected Filling
% =========================================================================

function [ebsd_pro, phaseMapClean, oriQuatClean] = doProtectedFilling(ebsd, oriQuatClean, phaseMapClean, protectedMask, sampleMask)  
% DOPROTECTEDFILLING Fills pixels previously flagged as protected (high MAD/WSR) using specific rules.
    global params
    radii = params.radius_fill; 
    coverageFrac = params.coverageFrac;
    thresholdFracRing = params.thresholdFracRing; 
    misTol_ori = params.misTol_ori; 
    % EBSD data handling
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol] = size(phaseMapClean);
    ebsdPhaseIdGrid = reshape(ebsd.phaseId, Nrow, Ncol);
    % Initial mask of pixels to be filled
    initialProtectedHoleMask = (ebsdPhaseIdGrid == notIndexedId) & protectedMask & sampleMask;
    holeMask = initialProtectedHoleMask;
    
    % Loop over radii in reverse
    fprintf('\n--- Starting Protected Hole Filling ---\n');
    fillCountTotal = 0;
    
    for k = numel(radii):-1:1
        radius = radii(k);
        fprintf('Processing Radius = %d ...\n', radius);
        kernelFull = double(fspecial('disk', radius)) > 0;
        N = size(kernelFull);
        innerKernel = false(N); 
        innerKernel(radius+1, radius+1) = true; 
        outerKernel = kernelFull & ~innerKernel;             
        fillCountRadius = 0;
        passes = 0;
        baseHoleMask = holeMask;
        % track fills per pass
        filledLastPass = false(Nrow, Ncol);
        % Multiple passes
        while true
            passes = passes + 1;           
            % --- Initialise padded maps for current pass ---
            paddedPhase = padarray(phaseMapClean, [radius radius], 0, 'both');
            paddedOri   = padarray(oriQuatClean, [radius radius], 0, 'both');
            
            % Neighbour-based selection
            if passes == 1
                currentHoleMask = baseHoleMask & (phaseMapClean == notIndexedId);
            else
                dilatedFilled = imdilate(filledLastPass, outerKernel);
                currentHoleMask = baseHoleMask & dilatedFilled & (phaseMapClean == notIndexedId);
            end
                       
            [holeRows, holeCols] = find(currentHoleMask);
            numHoles = numel(holeRows);
            
            if numHoles == 0
                fprintf('⚠ No remaining protected holes to fill at this radius.\n');
                break;
            end
            
            fillCountPass = 0;
            fprintf('Total candidate protected holes: %d\n', numHoles);
            
            % Loop over each hole
            for h = 1:numHoles
                i = holeRows(h); j = holeCols(h);
                if phaseMapClean(i,j) ~= notIndexedId
                    continue;
                end             
                iP = i + radius; jP = j + radius;               
                % Extract local window
                winPhase = paddedPhase(iP-radius:iP+radius, jP-radius:jP+radius);
                winOri   = paddedOri(iP-radius:iP+radius, jP-radius:jP+radius, :);                
                % Identify valid neighbours
                maskValid = (winPhase > 0 & winPhase ~= notIndexedId) & kernelFull;
                neighPhases = winPhase(maskValid);
                if isempty(neighPhases), continue; end
                
                % --- Check 1: Coverage ---
                coverage = numel(neighPhases) / sum(kernelFull(:));
                if coverage < coverageFrac
                    continue;
                end
                
                % --- Check 2: Phase Uniformity (all indexed neighbours must be the same phase) ---
                uniquePhases = unique(neighPhases);
                if numel(uniquePhases) > 1
                    continue;
                end
                domPhase = uniquePhases(1);
                % Build quaternion list
                qList = reshape(winOri, [], 4);
                qList = qList(maskValid(:), :); 
                
                % --- Check 3: Orientation (Ring-based) ---
                meanOri = calc_mean_ori_ring(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, maskValid, thresholdFracRing);
                % --- Skip if no orientation was returned ---
                if isempty(meanOri)
                    continue;
                end
                
                % ===== SUCCESSFUL FILL =====
                phaseMapClean(i,j) = domPhase;
                oriQuatClean(i,j,:) = [meanOri.a meanOri.b meanOri.c meanOri.d];
                
                fillCountPass  = fillCountPass + 1;
                fillCountRadius = fillCountRadius + 1;
                fillCountTotal  = fillCountTotal + 1;
                
                paddedPhase(iP,jP) = domPhase;
                paddedOri(iP,jP,:) = [meanOri.a meanOri.b meanOri.c meanOri.d]; 
                
                % Mark this pixel as newly filled
                filledLastPass(i,j) = true;
            end
            
            fprintf('  Pass %d complete: filled %d holes.\n', passes, fillCountPass);
            dilatedFilled = imdilate(filledLastPass, outerKernel);
            holeMask = baseHoleMask & dilatedFilled;
            if fillCountPass == 0
                break;
            end
        end % end while (MultiPass)
        fprintf('  ✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
    end % end loop over radii
    fprintf('\n--- Protected Hole Filling Complete. Total filled: %d ---\n', fillCountTotal);
    
    % --- Rebuild EBSD object with filled data ---
    ebsd_pro = ebsd;
    ebsd_pro.phaseId(:) = phaseMapClean(:); 
    
    qFull_fill = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (ebsd_pro.phaseId == pid) & initialProtectedHoleMask(:);
        
        if ~any(mask), continue; end
        ebsd_pro(mask).orientations = orientation(qFull_fill(mask), ebsd(pname).CS);
    end
    ebsd_fill = ebsd_pro; % Use the final variable name used in the main script
end

function meanOri = calc_mean_ori_ring(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, maskValid, thresholdFrac)
% CALC_MEAN_ORI_RING Calculates mean orientation based on the uniformity of the nearest orientation ring.
    meanOri = [];
    distMap = bwdist(innerKernel); % distance from center
    distMap = round(distMap);
    distMap(~outerKernel) = 0;
    found = false;    
    
    % Loop over rings from inner to outer
    for ringWidth = 1:radius
        ringMask = (distMap == ringWidth);
        if ~any(ringMask(:)), continue; end        
        % Find pixels that are valid and in this ring
        validLinearIdx = find(maskValid);
        ringLinearIdx  = find(ringMask(:));
        [~, commonIdx, ~] = intersect(validLinearIdx, ringLinearIdx);        
        
        % Filter to dominant phase in this ring
        qRing = qList(commonIdx, :);
        neighPhasesRing = neighPhases(commonIdx);
        qRingDom = qRing(neighPhasesRing == domPhase, :);
        Nring = size(qRingDom, 1);
        
        % Check minimum number of neighbours
        if Nring < 2
            found = false;  % insufficient neighbours
            break;
        end        
        
        % Cluster orientations
        D_ring = 2*acos(min(abs(qRingDom*qRingDom.'),1));
        D_ring(1:size(D_ring,1)+1:end) = 0;
        Dc_ring = squareform(D_ring);
        Z_ring = linkage(Dc_ring,'single');
        idx_ring = cluster(Z_ring,'cutoff',misTol_ori,'criterion','distance');
        counts_ring = accumarray(idx_ring,1);
        [maxCountRing, domClusterRing] = max(counts_ring);
        domFracRing = maxCountRing / Nring;
        
        % Check threshold
        if domFracRing < thresholdFrac
            found = false;  
            break;
        else
            found = true;    
        end
    end    
    
    % If the last ring checked passed the threshold, compute the mean
    if found
        % Select only quaternions in the dominant orientation cluster
        membersRing = (idx_ring == domClusterRing);   % pick largest cluster
        qClusterRing = qRingDom(membersRing, :);
        % Compute mean orientation using only dominant cluster
        qAllDom = quaternion(qClusterRing.');
        meanOri = mean(qAllDom, 'meanOrientation');
    end
end

%% =========================================================================
%   7. General Utilities (Stats, Plots, Residuals)
% =========================================================================

function showPhaseStats(ebsdObj, phases, msg)
% SHOWPHASESTATS Prints a summary of phase distribution and coverage to the console.
    fprintf('\n%s\n', msg); fprintf('--------------------------------\n');
    total = numel(ebsdObj);
    for i=1:numel(phases)
        n=numel(ebsdObj(phases{i}));
        fprintf('%-12s: %6d points (%.2f%%)\n', phases{i}, n, 100*n/total);
    end
    fprintf('--------------------------------\n');
end

function plotPhaseMap(ebsdObj, sampleName, exportPath, suffix,res)
% PLOTPHASEMAP Generates a phase map figure and exports it to PNG.
    f=figure('Visible','off'); plot(ebsdObj,'phase');
    leg=legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','on','FontSize',10);
    leg.Position(1)=0.5-leg.Position(3)/2;
    savePNG(f,sprintf('%s_PhaseMap_%s',sampleName,suffix),exportPath,res);
end

function plotIPFMapPhases(ebsdObj, sampleName, exportPath, suffix, res)
% PLOTIPFMAPPHASES Generates individual Inverse Pole Figure (IPF) maps for each indexed phase and exports them.
    phases = ebsdObj.mineralList;
    for i = 1:numel(phases)
        pname = phases{i}; if strcmpi(pname,'notIndexed'), continue; end
        ebsdPhase = ebsdObj(pname); if isempty(ebsdPhase), continue; end
        f = figure('Visible','off');
        plot(ebsdPhase, ebsdPhase.orientations); axis equal; 
        savePNG(f, sprintf('%s_IPFMap_%s_%s', sampleName, pname, suffix), exportPath, res);
    end
end

function savePNG(figHandle,filenameStem,exportPath,res)
% SAVEPNG Exports a figure handle to a PNG file at a specified resolution.
    exportgraphics(figHandle,fullfile(exportPath,[filenameStem,'.png']),'Resolution',res);
    close(figHandle);
    fprintf('Saved: %s.png\n', filenameStem);
end
