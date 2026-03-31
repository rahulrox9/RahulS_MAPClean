%% MAPClean – Microstructurally Adaptive Pixel-Level Cleaning
% MAPClean is a phase- and orientation-aware EBSD cleaning workflow for
% sparsely indexed and noisy maps. The pipeline consists of:
%   1. MAD-based removal of unreliable pixels
%   2. Sample mask generation
%   3. Phase-wise statistical reassignment (WSR)
%   4. Orientation-wise spike removal / reassignment (WSR)
%   5. Hole filling using either BFS (strict mode) or MPF (relaxed mode)
%
% The workflow automatically switches between strict and relaxed filling
% regimes depending on the notIndexed fraction within the sample mask.
%
% Main outputs:
%   - cleaned EBSD object
%   - phase map (phaseMapClean)
%   - quaternion orientation map (oriQuatClean)
%   - checkpointed intermediate results

clc; clear; close all;
addpath(genpath(pwd));
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

c = parcluster('local');
maxWorkers = c.NumWorkers;
saveProfile(c);
delete(c.Jobs);
delete(gcp('nocreate'));
parpool('local', maxWorkers);


%% User Configuration
%% --- Stage Control ---
runStart    = true;
runMAD      = true;
runCrop     = true;
runPhaseWSR = true;
runOriWSR   = true;
runHoleFill = true;
runSaveFile = true;

%% --- Parameters ---
global params
params.exportRes      = 300;
params.madThreshold   = 1.0;
params.qcfrac         = 0.6;
params.radius_phase   = 3;
params.min_dom_frac   = 0.5;
params.misTol_ori     = 5*degree;
params.thresholdFrac  = 0.75;
params.minFrac_ori    = 0.25;
params.radius_ori     = 2;
params.minLead        = 2;
params.scaleLead      = 0.1;

params.radius_fill_strict   = [6 5 4 3 2 1];
params.phaseFrac_strict     = containers.Map('KeyType','double','ValueType','any');
% Format: [Target_Count / N_disk, Threshold / Target_Count]
params.phaseFrac_strict(6) = [62/136 45/62];
params.phaseFrac_strict(5) = [45/100 34/45];
params.phaseFrac_strict(4) = [30/68 23/30];
params.phaseFrac_strict(3) = [22/44 17/22];
params.phaseFrac_strict(2) = [13/20 11/13];
params.phaseFrac_strict(1) = [6/8 6/6];

params.radius_fill_relaxed  = [7 6 5 4 3 2 1];
params.phaseFrac_relaxed    = containers.Map('KeyType','double','ValueType','any');
params.phaseFrac_relaxed(7) = [72/184 56/72];
params.phaseFrac_relaxed(6) = [50/136 36/50];
params.phaseFrac_relaxed(5) = [34/100 26/34];
params.phaseFrac_relaxed(4) = [24/68 18/24];
params.phaseFrac_relaxed(3) = [16/44 12/16];
params.phaseFrac_relaxed(2) = [10/20 8/10];
params.phaseFrac_relaxed(1) = [4/8 4/4];

%% Directories
dataDir = fullfile(pwd,'DataFiles');
sampleNames = {'sample1','sample2'}; # samplename must match the "filename.ctf" without the ".ctf"
fileList = cellfun(@(s) dir(fullfile(dataDir,[s '.ctf'])), sampleNames, 'UniformOutput', false);
fileList = vertcat(fileList{:});
checkpointDir = fullfile(pwd,'checkpoints');
exportDir = fullfile(pwd,'MAPCLean');
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end
if ~exist(exportDir,'dir'), mkdir(exportDir); end

%% Loop Samples
for fi = 1:numel(fileList)
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    exportPath = fullfile(exportDir, sampleName);
    if ~exist(exportPath,'dir'), mkdir(exportPath); end

    diaryFile = fullfile(exportDir, [sampleName '_logfile.txt']);
    diary off;
    if exist(diaryFile,'file'), delete(diaryFile); end
    diary(diaryFile); diary on;

    fprintf('\n===== Processing Sample: %s =====\n', sampleName);

    ebsd_raw = EBSD.load(fullfile(dataDir,fileList(fi).name), ...
        'convertSpatial2EulerReferenceFrame').gridify;

    phases        = ebsd_raw.mineralList;
    notIndexedId  = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds   = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol]  = ebsd_raw.size;

    ForsteriteColor = str2rgb('red');
    DiopsideColor   = str2rgb('blue');
    AnorthiteColor  = str2rgb('yellow');
    for i = 1:numel(MinPhaseNames)
        pname = MinPhaseNames{i};
        if strcmpi(pname,'Forsterite')
            ebsd_raw(pname).color = ForsteriteColor;
        elseif strcmpi(pname,'Diopside')
            ebsd_raw(pname).color = DiopsideColor;
        elseif strcmpi(pname,'Anorthite')
            ebsd_raw(pname).color = AnorthiteColor;
        end
    end
    clear ForsteriteColor DiopsideColor AnorthiteColor

    %% --- Checkpoint Paths ---
    madFile   = fullfile(checkpointDir, sprintf('%s_ebsd_mad.mat', sampleName));
    cropFile  = fullfile(checkpointDir, sprintf('%s_crop.mat', sampleName));
    phaseFile = fullfile(checkpointDir, sprintf('%s_ebsd_phase.mat', sampleName));
    oriFile   = fullfile(checkpointDir, sprintf('%s_ebsd_ori.mat', sampleName));
    fillFile  = fullfile(checkpointDir, sprintf('%s_ebsd_fill.mat', sampleName));
    paramFile = fullfile(checkpointDir, sprintf('%s_params.mat', sampleName));
    
    ebsd = ebsd_raw;
    % START
    if runStart
        tStartStats = tic;
        showPhaseStats(ebsd_raw, phases, 'Phase distribution before cleaning');
        fprintf('✔ Raw stats done (%.2f s)\n', toc(tStartStats));

        tStartPlot = tic;
        plotPhaseMap(ebsd_raw, sampleName, exportPath, 'raw', params.exportRes);
        plotIPFMapPhases(ebsd_raw, sampleName, exportPath, 'raw', params.exportRes);
        fprintf('✔ Raw figure export done (%.2f s)\n', toc(tStartPlot));
    end
    % MAD FILTER
    tMAD = tic;
    fprintf('--- MAD Filter ---\n');
    if runMAD
        [ebsd_mad, badPixels] = doMADFilter(ebsd, sampleName, exportPath);
        save(madFile,'ebsd_mad','badPixels');
        ebsd = ebsd_mad;
        fprintf('✔ MAD checkpoint saved (%.2f s)\n', toc(tMAD));
    elseif exist(madFile,'file')
        load(madFile,'ebsd_mad','badPixels');
        ebsd = ebsd_mad;
        fprintf('✔ MAD checkpoint loaded (%.2f s)\n', toc(tMAD));
    else
        badPixels = false(Nrow,Ncol);
        fprintf('⚠ MAD skipped and no checkpoint exists (%.2f s)\n', toc(tMAD));
    end
    % SAMPLE CROP
    tCrop = tic;
    fprintf('--- Sample Mask ---\n');
    if runCrop
        [sampleMask, indexedMask] = buildSampleMask(ebsd);
        save(cropFile,'sampleMask','indexedMask');
        fprintf('✔ Sample mask saved (%.2f s)\n', toc(tCrop));
        tCropStats = tic;
        showPhaseStats(ebsd(sampleMask), ebsd.mineralList, 'Phase distribution after cropping');
        fprintf('✔ Sample-mask stats done (%.2f s)\n', toc(tCropStats));
    elseif exist(cropFile,'file')
        load(cropFile,'sampleMask','indexedMask');
        fprintf('✔ Sample mask checkpoint loaded (%.2f s)\n', toc(tCrop));
    else
        sampleMask = true(Nrow, Ncol);
        indexedMask = sampleMask;
        fprintf('⚠ Crop skipped and no checkpoint exists (%.2f s)\n', toc(tCrop));
    end
    % EVALUATION - RELAXED OR STRICT? 
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    validPhaseIds = ebsd.phaseId(sampleMask);
    fracNotIndexed = sum(validPhaseIds == notIndexedId) / numel(validPhaseIds);
    if fracNotIndexed < params.qcfrac
        runstrict = true;
        params.radius_fill = params.radius_fill_strict;
        params.phaseFrac   = params.phaseFrac_strict;
        fprintf('✔ Data sufficiently indexed (%.2f%% notIndexed). strict mode\n', fracNotIndexed*100);
    else
        runstrict = false;
        params.radius_fill = params.radius_fill_relaxed;
        params.phaseFrac   = params.phaseFrac_relaxed;
        fprintf('⚠ Sparse data (%.2f%% notIndexed). relaxed mode\n', fracNotIndexed*100);
    end
    % WSR
    fprintf('--- WSR ---\n');
    if runstrict
        if runPhaseWSR
            [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR_strict(ebsd, badPixels, sampleName, exportPath, sampleMask);
            save(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_phase;
            fprintf('✔ strict phase WSR saved\n');
        elseif exist(phaseFile,'file')
            load(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_phase;
            fprintf('✔ strict phase WSR loaded\n');
        end

        if runOriWSR
            [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath, sampleMask);
            save(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_ori;
            fprintf('✔ Orientation WSR saved\n');
        elseif exist(oriFile,'file')
            load(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_ori;
            fprintf('✔ Orientation WSR loaded\n');
        end
    else
        if runPhaseWSR
            [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR_relaxed(ebsd, badPixels, sampleName, exportPath, sampleMask);
            save(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_phase;
            fprintf('✔ relaxed phase WSR saved\n');
        elseif exist(phaseFile,'file')
            load(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_phase;
            fprintf('✔ relaxed phase WSR loaded\n');
        end

        if runOriWSR
            [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath, sampleMask);
            save(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_ori;
            fprintf('✔ Orientation WSR saved\n');
        elseif exist(oriFile,'file')
            load(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean','protectedMask');
            ebsd = ebsd_ori;
            fprintf('✔ Orientation WSR loaded\n');
        end
    end
    % HOLE FILLING
    tFill = tic;
    fprintf('--- Filling Holes ---\n');
    if runHoleFill
        if runstrict
            [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingBFS(ebsd, oriQuatClean, phaseMapClean, protectedMask, sampleMask);
        else
            [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingMPF(ebsd, oriQuatClean, phaseMapClean, protectedMask, sampleMask);
        end
        ebsd = EBSD(ebsd_fill, 'convert');
        fprintf('✔ Hole filling core done (%.2f s)\n', toc(tFill));
    
        tFillStats = tic;
        showPhaseStats(ebsd, ebsd.mineralList, 'Phase distribution after hole fill');
        fprintf('✔ Hole-fill stats done (%.2f s)\n', toc(tFillStats));
    
        tFillPlot = tic;
        plotPhaseMap(ebsd, sampleName, exportPath, 'fill', params.exportRes);
        plotIPFMapPhases(ebsd, sampleName, exportPath, 'fill', params.exportRes);
        fprintf('✔ Hole-fill figure export done (%.2f s)\n', toc(tFillPlot));
    
        save(fillFile,'ebsd_fill','phaseMapClean','oriQuatClean','protectedMask');
        fprintf('✔ Hole fill checkpoint saved (%.2f s)\n', toc(tFill));
    elseif exist(fillFile,'file')
        load(fillFile,'ebsd_fill','phaseMapClean','oriQuatClean','protectedMask');
        ebsd = ebsd_fill;
        fprintf('✔ Hole fill checkpoint loaded (%.2f s)\n', toc(tFill));
    else
        fprintf('⚠ Hole fill skipped and no checkpoint exists (%.2f s)\n', toc(tFill));
    end
    % SAVE CLEAN FILE
    if runSaveFile
        outFile = fullfile(dataDir,[sampleName '_clean.ctf']);
        ebsd.export(outFile);
        save(paramFile,'params');
        fprintf('✔ Saved cleaned EBSD: %s\n', outFile);
        fprintf('✔ Parameters saved\n');
    end

    diary off;
end
%% =========================================================================
%  MAD FILTER
%% =========================================================================
function [ebsd_mad, badPixels] = doMADFilter(ebsd, sampleName, exportPath)
    % Purpose:
    %   Removes unreliable EBSD pixels using a threshold on the MAD value and
    %   reassigns them to notIndexed.
    %
    % Inputs:
    %   ebsd        - EBSD map
    %   sampleName  - sample identifier for export
    %   exportPath  - output directory for figures
    %
    % Outputs:
    %   ebsd_mad    - EBSD map after MAD filtering
    %   badPixels   - logical mask of removed pixels
    %
    % Logic:
    %   Pixels with MAD > params.madThreshold are considered unreliable and are
    %   reset to notIndexed. Phase statistics and phase/IPF plots are exported.
    global params
    tMAD = tic;
    fprintf('\n--- Applying MAD Filter (Threshold = %.2f rad) ---\n', params.madThreshold);

    badPixels = ebsd.mad > params.madThreshold;
    numBad = sum(badPixels,'all');
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    ebsd(badPixels).phaseId = notIndexedId;
    ebsd_mad = ebsd;
    fprintf('✔ MAD filter applied: %d pixels set to notIndexed (%.2f s)\n', numBad, toc(tMAD));

    tStats = tic;
    showPhaseStats(ebsd_mad, phases, 'Phase distribution after MAD filter');
    fprintf('✔ MAD stats done (%.2f s)\n', toc(tStats));

    tPlot = tic;
    plotPhaseMap(ebsd_mad, sampleName, exportPath, 'MADfilter', params.exportRes);
    plotIPFMapPhases(ebsd_mad, sampleName, exportPath, 'MADfilter', params.exportRes);
    fprintf('✔ MAD figure export done (%.2f s)\n', toc(tPlot));

    fprintf('✔ MAD total done (%.2f s)\n', toc(tMAD));
end
%% =========================================================================
%  SAMPLE MASK
%% =========================================================================
function [sampleMask, indexedMask] = buildSampleMask(ebsd)
    % Purpose:
    %   Constructs a sample support mask from indexed pixels so that later
    %   cleaning and filling operations are restricted to the mapped specimen.
    %
    % Inputs:
    %   ebsd - EBSD map
    %
    % Outputs:
    %   sampleMask  - logical mask covering the interior of the specimen
    %   indexedMask - logical mask of currently indexed pixels
    %
    % Logic:
    %   The sample mask is formed row-wise and column-wise by spanning between
    %   the first and last indexed pixel in each row/column.
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    [Nrow, Ncol] = ebsd.size;
    indexedMask = reshape(ebsd.phaseId ~= notIndexedId, Nrow, Ncol);
    sampleMask = false(Nrow, Ncol);
    for r = 1:Nrow
        idx = find(indexedMask(r,:));
        if isempty(idx), continue; end
        sampleMask(r, idx(1):idx(end)) = true;
    end
    for c = 1:Ncol
        idx = find(indexedMask(:,c));
        if isempty(idx), continue; end
        sampleMask(idx(1):idx(end), c) = true;
    end
end
%% =========================================================================
%  PHASE WSR — STRICT
%% =========================================================================
function [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR_strict(ebsd, badPixels, sampleName, exportPath, sampleMask)
    % Purpose:
    %   Performs strict phase-wise statistical reassignment (WSR) using local
    %   phase neighbourhood support and limited orientation updating.
    %
    % Inputs:
    %   ebsd        - EBSD map
    %   badPixels   - MAD-filtered pixel mask
    %   sampleName  - sample identifier
    %   exportPath  - figure export directory
    %   sampleMask  - valid specimen mask
    %
    % Outputs:
    %   ebsd_phase     - EBSD map after strict phase reassignment
    %   phaseMapClean  - cleaned phase map
    %   oriQuatClean   - updated quaternion orientation map
    %   protectedMask  - pixels protected from later ordinary hole filling
    %
    % Logic:
    %   Local phase support is evaluated within a radius-dependent disk
    %   neighbourhood. Weakly supported pixels may be removed, while pixels with
    %   sufficiently dominant neighbouring phase support may be reassigned.
    global params
    tWSR = tic;
    fprintf('Starting Phase WSR (strict, Radius = %d)\n', params.radius_phase);
    phases = ebsd.mineralList;
    Nrow = ebsd.size(1);
    Ncol = ebsd.size(2);
    phaseMapoG = reshape(double(ebsd.phaseId), Nrow, Ncol);
    phaseMapClean = phaseMapoG;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    oriQuatClean = zeros(Nrow,Ncol,4);
    N_total = Nrow*Ncol;

    % Protected mask = original MAD-removed pixels + later indecisive removals
    protectedMask = badPixels;

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

    kernel_phase = double(fspecial('disk', params.radius_phase)) > 0;
    kernel_phase(params.radius_phase+1, params.radius_phase+1) = 0;

    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        fprintf('Processing phase: %s\n', pname);
        mask = (phaseMapoG == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        [rows, cols] = find(mask);
        paddedPhase = padarray(phaseMapClean,[params.radius_phase params.radius_phase],0,'both');
        paddedOri   = padarray(oriQuatClean,[params.radius_phase params.radius_phase],0,'both');

        for k = 1:numel(rows)
            i = rows(k);
            j = cols(k);
            iP = i + params.radius_phase;
            jP = j + params.radius_phase;
            win = paddedPhase(iP-params.radius_phase:iP+params.radius_phase, ...
                              jP-params.radius_phase:jP+params.radius_phase);
            neigh = win(kernel_phase);
            total_neigh = numel(neigh);
            indexedMask = (neigh > 0 & neigh ~= notIndexedId);
            validNeigh = neigh(indexedMask);
            Ni = numel(validNeigh) / total_neigh;

            if Ni <= 0.25
                phaseMapClean(i,j) = notIndexedId;
                oriQuatClean(i,j,:) = 0;
                protectedMask(i,j) = true;
                continue;
            end

            if Ni > 0.25 && Ni < 0.5
                uniquePhases = unique(validNeigh);
                if isscalar(uniquePhases)
                    maj = uniquePhases;
                    if phaseMapClean(i,j) ~= maj
                        phaseMapClean(i,j) = maj;
                    end
                    oriWin = paddedOri(iP-params.radius_phase:iP+params.radius_phase, ...
                                       jP-params.radius_phase:jP+params.radius_phase,:);
                    oriList = reshape(oriWin, [], 4);
                    validOriMask = (win(kernel_phase) == maj);
                    neighbourQuatsList = oriList(validOriMask,:);
                    Nneighbours = size(neighbourQuatsList,1);
                    if Nneighbours >= 2
                        currentQ_vec = squeeze(oriQuatClean(i,j,:))';
                        [meanOri, ~] = calc_mean_ori_wsr_strict(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                        oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                    end
                end
                continue;
            end

            fracThresh = (Ni >= 0.75) * (2/3) + (Ni < 0.75) * (3/4);
            maj = mode(validNeigh);
            fracMaj = sum(validNeigh == maj) / numel(validNeigh);
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
                    [meanOri, ~] = calc_mean_ori_wsr_strict(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                    oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                end
            end
        end
    end

    ebsd_phase = ebsd;
    ebsd_phase.phaseId(:) = phaseMapClean(:);
    qFull_phase = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        ebsd_phase(mask).orientations = orientation(qFull_phase(mask), ebsd(pname).CS);
    end
    fprintf('✔ strict Phase WSR core done (%.2f s)\n', toc(tWSR));
    tStats = tic;
    showPhaseStats(ebsd_phase, phases, 'Phase distribution after Phase WSR');
    fprintf('✔ strict Phase WSR stats done (%.2f s)\n', toc(tStats));
    tPlot = tic;
    plotPhaseMap(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    plotIPFMapPhases(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    fprintf('✔ strict Phase WSR figure export done (%.2f s)\n', toc(tPlot));
    fprintf('✔ strict Phase WSR total done (%.2f s)\n', toc(tWSR));
end
%% =========================================================================
%  PHASE WSR — RELAXED
%% =========================================================================
function [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR_relaxed(ebsd, badPixels, sampleName, exportPath, sampleMask)
    % Purpose:
    %   Performs relaxed phase-wise statistical reassignment (WSR) for sparse
    %   datasets where strict support requirements would remove too many pixels.
    %
    % Inputs:
    %   ebsd        - EBSD map
    %   badPixels   - MAD-filtered pixel mask
    %   sampleName  - sample identifier
    %   exportPath  - figure export directory
    %   sampleMask  - valid specimen mask
    %
    % Outputs:
    %   ebsd_phase     - EBSD map after strict phase reassignment
    %   phaseMapClean  - cleaned phase map
    %   oriQuatClean   - updated quaternion orientation map
    %   protectedMask  - pixels protected from later ordinary hole filling
    %
    % Logic:
    %   Pixels with extremely poor support are removed, while pixels whose local
    %   neighbourhood shows sufficient dominant-phase agreement may be flipped
    %   to that phase using a lower support threshold than in strict mode.
    global params
    tWSR = tic;
    minDomFrac = params.min_dom_frac;
    fprintf('\n--- Starting Phase WSR (relaxed, Radius = %d) ---\n', params.radius_phase);
    phases = ebsd.mineralList;
    Nrow = ebsd.size(1);
    Ncol = ebsd.size(2);
    phaseMapoG = reshape(double(ebsd.phaseId), Nrow, Ncol);
    phaseMapClean = phaseMapoG;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    oriQuatClean = zeros(Nrow,Ncol,4);
    N_total = Nrow*Ncol;

    protectedMask = badPixels;

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

    kernel_phase = double(fspecial('disk', params.radius_phase)) > 0;
    kernel_phase(params.radius_phase+1, params.radius_phase+1) = 0;

    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        fprintf('Processing phase: %s\n', pname);
        mask = (phaseMapoG == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        [rows, cols] = find(mask);
        paddedPhase = padarray(phaseMapClean,[params.radius_phase params.radius_phase],0,'both');
        paddedOri   = padarray(oriQuatClean,[params.radius_phase params.radius_phase],0,'both');
        countRemoved = 0;
        countPhaseFlip = 0;

        for k = 1:numel(rows)
            i = rows(k);
            j = cols(k);
            iP = i + params.radius_phase;
            jP = j + params.radius_phase;
            win = paddedPhase(iP-params.radius_phase:iP+params.radius_phase, ...
                              jP-params.radius_phase:jP+params.radius_phase);
            neigh = win(kernel_phase);
            total_neigh = numel(neigh);
            indexedMask = (neigh > 0 & neigh ~= notIndexedId);
            validNeigh = neigh(indexedMask);
            Ni = numel(validNeigh) / total_neigh;

            if Ni <= 0.10
                phaseMapClean(i,j) = notIndexedId;
                oriQuatClean(i,j,:) = 0;
                protectedMask(i,j) = true;
                countRemoved = countRemoved + 1;
                continue;
            end

            [uniquePh, ~, ic] = unique(validNeigh);
            counts = accumarray(ic,1);
            [maxCount, idxMax] = max(counts);
            domPhase = uniquePh(idxMax);
            fracMaj = maxCount / numel(validNeigh);
            if phaseMapClean(i,j) ~= domPhase && fracMaj >= minDomFrac
                phaseMapClean(i,j) = domPhase;
                countPhaseFlip = countPhaseFlip + 1;
                oriWin = paddedOri(iP-params.radius_phase:iP+params.radius_phase, ...
                                   jP-params.radius_phase:jP+params.radius_phase,:);
                oriList = reshape(oriWin, [], 4);
                validOriMask = (win(kernel_phase) == domPhase);
                neighbourQuatsList = oriList(validOriMask,:);
                Nneighbours = size(neighbourQuatsList,1);
                if Nneighbours >= 2
                    currentQ_vec = squeeze(oriQuatClean(i,j,:))';
                    [meanOri, ~] = calc_mean_ori_wsr_relaxed(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                    oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                end
            end
        end

        fprintf('  ✔ Phase %s: Removed %d points. Flipped %d points.\n', pname, countRemoved, countPhaseFlip);
    end

    ebsd_phase = ebsd;
    ebsd_phase.phaseId(:) = phaseMapClean(:);
    qFull_phase = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid) & sampleMask;
        if ~any(mask,'all'), continue; end
        ebsd_phase(mask).orientations = orientation(qFull_phase(mask), ebsd(pname).CS);
    end
    fprintf('✔ relaxed Phase WSR core done (%.2f s)\n', toc(tWSR));
    tStats = tic;
    showPhaseStats(ebsd_phase, phases, 'Phase distribution after Phase WSR');
    fprintf('✔ relaxed Phase WSR stats done (%.2f s)\n', toc(tStats));
    tPlot = tic;
    plotPhaseMap(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    plotIPFMapPhases(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    fprintf('✔ relaxed Phase WSR figure export done (%.2f s)\n', toc(tPlot));
    fprintf('✔ relaxed Phase WSR total done (%.2f s)\n', toc(tWSR));
end
%% =========================================================================
%  ORIENTATION WSR 
%% =========================================================================
function [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath, sampleMask)
    % Purpose:
    %   Detects and corrects orientation outliers ("wild spikes") within each
    %   phase after phase cleaning.
    %
    % Inputs:
    %   ebsd          - EBSD map
    %   oriQuatClean  - quaternion orientation map
    %   phaseMapClean - cleaned phase map
    %   MinPhaseIds   - indexed phase IDs
    %   sampleName    - sample identifier
    %   exportPath    - figure export directory
    %   sampleMask    - valid specimen mask
    %
    % Outputs:
    %   ebsd_ori      - EBSD map after orientation WSR
    %   oriQuatClean  - updated quaternion map
    %   phaseMapClean - phase map (may include removed spikes)
    %
    % Logic:
    %   For each phase, pixels with insufficient local orientation similarity
    %   are flagged as spikes. Their orientations are reassigned using local
    %   clustering-based mean orientation estimates.
    global params
    t0 = tic;
    fprintf('\nOrientation WSR:\n');

    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseNames = phases(MinPhaseIds);

    r = params.radius_ori;
    kernel_ori = double(fspecial('disk', r)) > 0;
    kernel_ori(r+1, r+1) = 0;

    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};

        mask = (phaseMapClean == pid) & sampleMask;
        [rowsAll, colsAll] = find(mask);
        numPx = numel(rowsAll);

        paddedOri = padarray(oriQuatClean, [r r], 0, 'both');
        paddedPhase = padarray(phaseMapClean, [r r], 0, 'both');

        fprintf('  Phase %s: %d pixels — detecting wild spikes (parfor)...\n', pname, numPx);
        tSpike = tic;

        % --------------------------------------------------------------
        % Precompute neighbour quaternion lists outside parfor
        % --------------------------------------------------------------
        neighQuatLists = cell(numPx,1);
        currentQList = zeros(numPx,4);

        for k = 1:numPx
            i = rowsAll(k);
            j = colsAll(k);
            iP = i + r;
            jP = j + r;

            oriWin = paddedOri(iP-r:iP+r, jP-r:jP+r, :);
            phaseWin = paddedPhase(iP-r:iP+r, jP-r:jP+r);

            nMask = (phaseWin == pid) & kernel_ori;
            oriList = reshape(oriWin, [], 4);
            neighQuatLists{k} = oriList(nMask(:), :);

            currentQList(k,:) = squeeze(oriQuatClean(i,j,:))';
        end

        wildSpikeFlags = false(numPx, 1);
        misTol_local = params.misTol_ori;
        minFrac_ori_local = params.minFrac_ori;

        parfor k = 1:numPx
            nqList = neighQuatLists{k};
            if isempty(nqList)
                continue;
            end

            cQ = currentQList(k,:);
            dots = abs(nqList * cQ');
            dots(dots > 1) = 1;

            fracSim = sum(2*acos(dots) < misTol_local) / numel(dots);
            if fracSim < minFrac_ori_local
                wildSpikeFlags(k) = true;
            end
        end

        spikeIdx = find(wildSpikeFlags);
        numSpikes = numel(spikeIdx);
        fprintf('  Phase %s: %d wild spikes found (%.1f sec)\n', pname, numSpikes, toc(tSpike));

        numRemoved = 0;
        numReoriented = 0;
        numSkippedTwin = 0;

        for ks = 1:numSpikes
            k = spikeIdx(ks);
            i = rowsAll(k);
            j = colsAll(k);
            iP = i + r;
            jP = j + r;

            oriWin = paddedOri(iP-r:iP+r, jP-r:jP+r, :);
            phaseWin = paddedPhase(iP-r:iP+r, jP-r:jP+r);

            nMask = (phaseWin == pid) & kernel_ori;
            oriList = reshape(oriWin, [], 4);
            nqList = oriList(nMask(:), :);
            Nn = size(nqList, 1);

            cQ_vec = squeeze(oriQuatClean(i,j,:))';
            currentQ = quaternion(cQ_vec(1), cQ_vec(2), cQ_vec(3), cQ_vec(4));

            if Nn == 0
                phaseMapClean(i,j) = notIndexedId;
                oriQuatClean(i,j,:) = 0;
                paddedOri(iP,jP,:) = 0;
                numRemoved = numRemoved + 1;
                continue;
            end

            if Nn == 1
                continue;
            end

            if Nn == 2
                dots = abs(nqList * cQ_vec');
                dots(dots > 1) = 1;
                if all(2*acos(dots) < params.misTol_ori)
                    qM = mean(quaternion(nqList(:,1), nqList(:,2), nqList(:,3), nqList(:,4)), 'meanOrientation');
                    oriQuatClean(i,j,:) = [qM.a qM.b qM.c qM.d];
                    paddedOri(iP,jP,:) = [qM.a qM.b qM.c qM.d];
                    numReoriented = numReoriented + 1;
                end
                continue;
            end

            dots = abs(nqList * cQ_vec');
            dots(dots > 1) = 1;

            if all(2*acos(dots) < params.misTol_ori)
                qM = mean(quaternion(nqList.'), 'meanOrientation');
            else
                qM = calc_mean_ori_wsr_strict(nqList, params.misTol_ori, Nn, cQ_vec);
            end

            if strcmpi(pname,'Anorthite')
                cs = ebsd('Anorthite').CS;
                misOri = qM * inv(currentQ);
                twinAngleTol = 5 * degree;
                axisTol = 5 * degree;
                misAng = angle(misOri);
                isTwin = false;

                if abs(misAng - 180*degree) <= twinAngleTol
                    twinAxis = axis(misOri);

                    candidateDirs = [ ...
                        vector3d(Miller(0,1,0,cs)), ...
                        vector3d(Miller(0,1,0,cs,'uvw')), ...
                        vector3d(Miller(0,0,1,cs,'uvw')), ...
                        vector3d(Miller(0,0,1,cs)), ...
                        vector3d(Miller(0,2,1,cs)), ...
                        vector3d(Miller(0,-2,1,cs)), ...
                        vector3d(Miller(1,0,0,cs,'uvw')), ...
                        vector3d(Miller(1,0,0,cs)) ...
                    ];

                    for tt = 1:numel(candidateDirs)
                        d = candidateDirs(tt);
                        if angle(twinAxis, d) < axisTol || angle(twinAxis, -d) < axisTol
                            isTwin = true;
                            break;
                        end
                    end
                end

                if isTwin
                    numSkippedTwin = numSkippedTwin + 1;
                    continue;
                end
            end

            oriQuatClean(i,j,:) = [qM.a qM.b qM.c qM.d];
            paddedOri(iP,jP,:) = [qM.a qM.b qM.c qM.d];
            numReoriented = numReoriented + 1;
        end

        fprintf('  %s OriWSR: Removed=%d  Reoriented=%d  TwinsSkipped=%d\n', ...
            pname, numRemoved, numReoriented, numSkippedTwin);

        clear neighQuatLists currentQList paddedOri paddedPhase
    end

    ebsd_ori = ebsd;
    ebsd_ori.phaseId(:) = phaseMapClean(:);

    qFull = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid) & sampleMask;
        if ~any(mask,'all')
            continue;
        end
        ebsd_ori(mask).orientations = orientation(qFull(mask), ebsd_ori(pname).CS);
    end

    fprintf('✔ Orientation WSR core done (%.2f s)\n', toc(t0));

    tStats = tic;
    showPhaseStats(ebsd_ori, phases, 'After Orientation WSR');
    fprintf('✔ Orientation WSR stats done (%.2f s)\n', toc(tStats));

    tPlot = tic;
    plotPhaseMap(ebsd_ori, sampleName, exportPath, 'OriWSR', params.exportRes);
    plotIPFMapPhases(ebsd_ori, sampleName, exportPath, 'OriWSR', params.exportRes);
    fprintf('✔ Orientation WSR figure export done (%.2f s)\n', toc(tPlot));

    fprintf('✔ Orientation WSR total done (%.2f s)\n', toc(t0));
end
%% =========================================================================
%  CALC MEAN ORI WSR STRICT
%% =========================================================================
function [meanOri, clusterSizes] = calc_mean_ori_wsr_strict(qList, misTol_ori, Nneighbours, currentQ_vec)
    % Purpose:
    %   Computes a representative local mean orientation under strict clustering rules.
    %
    % Inputs:
    %   qList         - neighbouring quaternions
    %   misTol_ori    - misorientation cutoff
    %   Nneighbours   - number of neighbours
    %   currentQ_vec  - current pixel quaternion
    %
    % Outputs:
    %   meanOri       - selected mean orientation
    %   clusterSizes  - cluster population sizes
    %
    % Logic:
    %   Orientations are clustered by single-linkage misorientation distance.
    %   If one cluster is sufficiently dominant, its mean is used; otherwise the
    %   cluster closest to the current orientation is selected.
    global params
    thresholdFrac = params.thresholdFrac;
    currentQ = quaternion(currentQ_vec(1), currentQ_vec(2), currentQ_vec(3), currentQ_vec(4));
    D = 2 * acos(min(abs(qList*qList.'),1));
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);
    Z = linkage(Dc,'single');
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;
    domClusterFrac = maxCount / Nneighbours;
    if domClusterFrac >= thresholdFrac
        qCluster = quaternion(qList(members,:).');
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    end
    uniqueClusters = unique(idx);
    minMisorientation = inf;
    closestClusterMean = quaternion(0,0,0,1);
    for c = 1:numel(uniqueClusters)
        clusterMembers = (idx == uniqueClusters(c));
        qCluster = quaternion(qList(clusterMembers,:).');
        clusterMean = mean(qCluster, 'meanOrientation');
        mis = angle(clusterMean * currentQ);
        if mis < minMisorientation
            minMisorientation = mis;
            closestClusterMean = clusterMean;
        end
    end
    meanOri = closestClusterMean;
end
%% =========================================================================
%  CALC MEAN ORI WSR RELAXED
%% =========================================================================
function [meanOri, clusterSizes] = calc_mean_ori_wsr_relaxed(qList, misTol_ori, Nneighbours, currentQ_vec)
    % Purpose:
    %   Computes a representative local mean orientation under relaxed
    %   clustering rules for sparse or noisy support.
    %
    % Inputs:
    %   qList         - neighbouring quaternions
    %   misTol_ori    - misorientation cutoff
    %   Nneighbours   - number of neighbours
    %   currentQ_vec  - current pixel quaternion
    %
    % Outputs:
    %   meanOri       - selected mean orientation
    %   clusterSizes  - cluster population sizes
    %
    % Logic:
    %   Similar to the strict version, but the dominant cluster is accepted
    %   using a lead-based criterion rather than a fixed dominance fraction.
    global params
    minLead = params.minLead;
    scaleLead = params.scaleLead;
    currentQ = quaternion(currentQ_vec(1), currentQ_vec(2), currentQ_vec(3), currentQ_vec(4));
    D = 2 * acos(min(abs(qList*qList.'),1));
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);
    Z = linkage(Dc,'single');
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;
    otherCounts = counts;
    otherCounts(domCluster) = 0;
    leadDiff = maxCount - max(otherCounts, 0);
    reqLead = max(minLead, ceil(scaleLead * Nneighbours));
    if leadDiff >= reqLead
        qCluster = quaternion(qList(members,:).');
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    end
    uniqueClusters = unique(idx);
    minMisorientation = inf;
    closestClusterMean = quaternion(0,0,0,1);
    for c = 1:numel(uniqueClusters)
        clusterMembers = (idx == uniqueClusters(c));
        qCluster = quaternion(qList(clusterMembers,:).');
        clusterMean = mean(qCluster, 'meanOrientation');
        mis = angle(clusterMean * currentQ);
        if mis < minMisorientation
            minMisorientation = mis;
            closestClusterMean = clusterMean;
        end
    end
    meanOri = closestClusterMean;
end
%% =========================================================================
%  HOLE FILLING - STRICT (BFS)
%% =========================================================================
function [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingBFS(ebsd, oriQuatClean, phaseMapClean, protectedMask, sampleMask)
    % Purpose:
    %   Fills unprotected notIndexed holes using a parallel BFS-style frontier expansion strategy.
    %
    % Inputs:
    %   ebsd          - EBSD map
    %   oriQuatClean  - quaternion orientation map
    %   phaseMapClean - cleaned phase map
    %   protectedMask - pixels excluded from ordinary filling
    %   sampleMask    - valid specimen mask
    %
    % Outputs:
    %   ebsd_fill     - EBSD map after BFS filling
    %   phaseMapClean - updated phase map
    %   oriQuatClean  - updated quaternion map
    %
    % Logic:
    %   Holes are decomposed into connected components. Components that satisfy
    %   minimum neighbourhood support are scheduled to workers and filled by
    %   iterative frontier propagation. Phase acceptance is controlled by
    %   radius-specific support thresholds, and orientations are assigned using
    %   local clustering-based mean estimates.
    global params
    t0 = tic;
    fprintf('\nHole Filling BFS:\n');

    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases, 'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);

    [Nrow, Ncol] = size(phaseMapClean);
    misTol = params.misTol_ori;
    threshFrac = params.thresholdFrac;
    radii = params.radius_fill;

    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local');
    end
    maxWorkers = pool.NumWorkers;

    for radius = radii
        tRad = tic;
        fprintf('--- Hole Filling BFS (radius = %d) ---\n', radius);

        phaseFrac = params.phaseFrac(radius);
        kernelFull = double(fspecial('disk', radius)) > 0;
        outerK = kernelFull;
        outerK(radius+1, radius+1) = false;
        nTotal = sum(outerK(:));

        totalFilledRadius = 0;
        iteration = 0;

        while true
            iteration = iteration + 1;

            holeMask = (phaseMapClean == notIndexedId) & sampleMask & ~protectedMask;
            if ~any(holeMask(:))
                break;
            end

            cc = bwconncomp(holeMask, 8);
            numClusters = cc.NumObjects;
            clusterPixels = cc.PixelIdxList;
            clusterSizes = cellfun(@numel, clusterPixels);
            totalHolePx = sum(clusterSizes);

            minClusterSize = max(2, ceil(phaseFrac(1) * nTotal));
            validIdx = find(clusterSizes >= minClusterSize);

            if isempty(validIdx)
                fprintf('  Radius %d iter %d: 0/%d clusters (%d px) with valid Ni\n', ...
                    radius, iteration, numClusters, totalHolePx);
                break;
            end

            [~, sortOrder] = sort(clusterSizes(validIdx), 'descend');
            validIdx = validIdx(sortOrder);

            indexedBinary = double((phaseMapClean ~= notIndexedId) & sampleMask & ~protectedMask);
            neighbourCount = conv2(indexedBinary, double(outerK), 'same');
            Ni_map = neighbourCount / nTotal;

            niPass = false(numel(validIdx), 1);
            for k = 1:numel(validIdx)
                ci = validIdx(k);
                if any(Ni_map(clusterPixels{ci}) >= phaseFrac(1))
                    niPass(k) = true;
                end
            end
            validIdx = validIdx(niPass);
            numValid = numel(validIdx);

            fprintf('  Radius %d iter %d: %d/%d clusters (%d px) with valid Ni\n', ...
                radius, iteration, numValid, numClusters, totalHolePx);

            if numValid == 0
                break;
            end

            clusterState = zeros(1, numClusters, 'uint8');
            clusterState(validIdx) = 1;

            fillableQueue = validIdx(:).';
            frozenList = [];
            frozenBlockedBy = false(numClusters, maxWorkers);

            futures = cell(maxWorkers, 1);
            workerFree = true(maxWorkers, 1);
            activeHaloMask = cell(maxWorkers, 1);
            workerClusterSize = zeros(maxWorkers, 1);
            workerStartTime = zeros(maxWorkers, 1);
            workerClusterId = zeros(maxWorkers, 1);

            freeWorkerQueue = 1:maxWorkers;

            numActive = 0;
            numAssigned = 0;
            cumulativeFilled = 0;
            numFilledIter = 0;

            while (numAssigned < numValid) || (numActive > 0)
                while ~isempty(freeWorkerQueue) && ~isempty(fillableQueue)
                    currentFree = freeWorkerQueue;
                    freeWorkerQueue = [];

                    assignedSomething = false;

                    for q = 1:numel(currentFree)
                        workerIdx = currentFree(q);
                        foundJob = false;

                        pos = 1;
                        while pos <= numel(fillableQueue)
                            ci = fillableQueue(pos);

                            if clusterState(ci) ~= 1
                                fillableQueue(pos) = [];
                                continue;
                            end

                            [isConflict, blockerVec] = getClusterBlockers(ci);

                            if isConflict
                                clusterState(ci) = 2;
                                frozenBlockedBy(ci, :) = blockerVec;
                                frozenList(end+1) = ci; 
                                fillableQueue(pos) = [];
                                continue;
                            end

                            pixList = clusterPixels{ci};
                            [rows_tmp, cols_tmp] = ind2sub([Nrow Ncol], pixList);

                            rMin_c = max(min(rows_tmp) - radius, 1);
                            rMax_c = min(max(rows_tmp) + radius, Nrow);
                            cMin_c = max(min(cols_tmp) - radius, 1);
                            cMax_c = min(max(cols_tmp) + radius, Ncol);

                            pPh = padarray(phaseMapClean(rMin_c:rMax_c, cMin_c:cMax_c), [radius radius], 0, 'both');
                            pOri = padarray(oriQuatClean(rMin_c:rMax_c, cMin_c:cMax_c, :), [radius radius], 0, 'both');
                            pProt = padarray(protectedMask(rMin_c:rMax_c, cMin_c:cMax_c), [radius radius], true, 'both');
                            pSample = padarray(sampleMask(rMin_c:rMax_c, cMin_c:cMax_c), [radius radius], false, 'both');

                            futures{workerIdx} = parfeval(pool, @fillClusterBFS_patchWorker, 1, ...
                                pixList, pPh, pOri, pProt, pSample, rMin_c, cMin_c, ...
                                Nrow, Ncol, radius, kernelFull, outerK, nTotal, ...
                                notIndexedId, phaseFrac, misTol, threshFrac);

                            hMask = false(Nrow, Ncol);
                            hMask(pixList) = true;
                            activeHaloMask{workerIdx} = imdilate(hMask, kernelFull);

                            workerFree(workerIdx) = false;
                            workerClusterSize(workerIdx) = numel(pixList);
                            workerStartTime(workerIdx) = now;
                            workerClusterId(workerIdx) = ci;

                            fillableQueue(pos) = [];
                            clusterState(ci) = 3;

                            numActive = numActive + 1;
                            numAssigned = numAssigned + 1;
                            assignedSomething = true;
                            foundJob = true;

                            fprintf('    assigned %d/%d -> worker %d | cluster size = %d | active = %d\n', ...
                                numAssigned, numValid, workerIdx, numel(pixList), numActive);
                            break;
                        end

                        if ~foundJob
                            freeWorkerQueue(end+1) = workerIdx;
                        end
                    end

                    if ~assignedSomething
                        break;
                    end
                end

                if numActive == 0
                    break;
                end

                busyWorkers = find(~workerFree);
                activeFutures = [futures{busyWorkers}];
                [finPos, result] = fetchNext(activeFutures);
                workerIdx = busyWorkers(finPos);
                ciDone = workerClusterId(workerIdx);

                lFilled = result.localFilled;
                nFilled = sum(lFilled);

                clusterTotal = workerClusterSize(workerIdx);
                clusterPct = 100 * nFilled / max(1, clusterTotal);
                workerTime = (now - workerStartTime(workerIdx)) * 86400;

                if nFilled > 0
                    filledPix = result.pixList(lFilled);
                    lPh = result.localPhase(lFilled);
                    lOri = result.localOri(lFilled, :);

                    for ki = 1:numel(filledPix)
                        [ii, jj] = ind2sub([Nrow Ncol], filledPix(ki));
                        phaseMapClean(ii, jj) = lPh(ki);
                        oriQuatClean(ii, jj, :) = reshape(lOri(ki, :), 1, 1, 4);
                    end
                end

                numFilledIter = numFilledIter + nFilled;
                cumulativeFilled = cumulativeFilled + nFilled;

                workerFree(workerIdx) = true;
                futures{workerIdx} = [];
                activeHaloMask{workerIdx} = [];
                workerClusterSize(workerIdx) = 0;
                workerStartTime(workerIdx) = 0;
                workerClusterId(workerIdx) = 0;
                numActive = numActive - 1;

                if ciDone > 0
                    clusterState(ciDone) = 4;
                end

                fprintf('    W%d | cluster = %d/%d (%.2f%%) | total = %d | active = %d | time = %.2f s\n', ...
                    workerIdx, nFilled, clusterTotal, clusterPct, cumulativeFilled, numActive, workerTime);

                freeWorkerQueue(end+1) = workerIdx;

                candIdx = find(clusterState == 2 & frozenBlockedBy(:, workerIdx).').';

                for kk = 1:numel(candIdx)
                    ci = candIdx(kk);

                    frozenBlockedBy(ci, workerIdx) = false;
                    [stillBlocked, blockerVec] = getClusterBlockers(ci);

                    if stillBlocked
                        frozenBlockedBy(ci, :) = blockerVec;
                    else
                        frozenBlockedBy(ci, :) = false;
                        clusterState(ci) = 1;
                        fillableQueue(end+1) = ci; 
                    end
                end

                if ~isempty(frozenList)
                    frozenList = frozenList(clusterState(frozenList) == 2);
                end

                activeIds = workerClusterId(workerClusterId > 0).';
                assert(isempty(intersect(fillableQueue, frozenList)), ...
                    'Bookkeeping error: fillable and frozen overlap.');
                assert(isempty(intersect(fillableQueue, activeIds)), ...
                    'Bookkeeping error: fillable and active overlap.');
                assert(isempty(intersect(frozenList, activeIds)), ...
                    'Bookkeeping error: frozen and active overlap.');
            end

            fprintf('  Radius %d iter %d: filled %d pixels\n', radius, iteration, numFilledIter);
            totalFilledRadius = totalFilledRadius + numFilledIter;

            if numFilledIter == 0
                fprintf('  Radius %d: converged after %d iterations\n', radius, iteration);
                break;
            end
        end

        fprintf('  Radius %d total: %d pixels filled (%.1f sec)\n', radius, totalFilledRadius, toc(tRad));
    end

    ebsd_fill = ebsd;
    ebsd_fill.phaseId(:) = phaseMapClean(:);

    qFull_fill = quaternion( ...
        oriQuatClean(:,:,1), ...
        oriQuatClean(:,:,2), ...
        oriQuatClean(:,:,3), ...
        oriQuatClean(:,:,4));

    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid) & sampleMask;
        if ~any(mask, 'all')
            continue;
        end
        ebsd_fill(mask).orientations = orientation(qFull_fill(mask), ebsd(pname).CS);
    end

    fprintf('✔ BFS complete (%.1f sec total)\n', toc(t0));

    function [tf, blockerVec] = getClusterBlockers(ci)
        blockerVec = false(1, maxWorkers);
        pix = clusterPixels{ci};

        for w = 1:maxWorkers
            if workerFree(w)
                continue;
            end
            lockedW = activeHaloMask{w};
            if ~isempty(lockedW) && any(lockedW(pix))
                blockerVec(w) = true;
            end
        end

        tf = any(blockerVec);
    end
end

function result = fillClusterBFS_patchWorker(pixList, patchPhase, patchOri, patchProtected, patchSample, rMin, cMin, Nrow, Ncol, r, kernelFull, outerK, nTotal, notIndexedId, phaseFrac, misTol, threshFrac)
    % Purpose:
    %   Worker routine for filling a single hole component during BFS filling.
    %
    % Logic:
    %   The component is processed iteratively in patch-local coordinates. A
    %   candidate pixel is filled only if it satisfies neighbourhood support,
    %   dominant-phase, and orientation-coherence criteria.
    nPix = numel(pixList);
    [rows_c, cols_c] = ind2sub([Nrow Ncol], pixList);

    rowOff = rMin - 1 - r;
    colOff = cMin - 1 - r;

    iL_vec = rows_c - rowOff;
    jL_vec = cols_c - colOff;

    N = 2*r + 1;
    innerK = false(N);
    innerK(r+1, r+1) = true;

    localPhase = zeros(nPix, 1, 'double');
    localOri = zeros(nPix, 4);
    localFilled = false(nPix, 1);
    skipMask = false(nPix, 1);

    madeChange = true;
    while madeChange
        madeChange = false;

        for ki = 1:nPix
            if localFilled(ki) || skipMask(ki)
                continue;
            end

            iL = iL_vec(ki);
            jL = jL_vec(ki);

            if patchProtected(iL, jL) || ~patchSample(iL, jL)
                skipMask(ki) = true;
                continue;
            end

            winPhase = patchPhase(iL-r:iL+r, jL-r:jL+r);
            winOri = patchOri(iL-r:iL+r, jL-r:jL+r, :);
            winProt = patchProtected(iL-r:iL+r, jL-r:jL+r);
            winSample = patchSample(iL-r:iL+r, jL-r:jL+r);

            validOuter = outerK & winSample & ~winProt;
            neighPh = winPhase(validOuter);
            validPh = neighPh(neighPh > 0 & neighPh ~= notIndexedId);
            Ni = numel(validPh) / nTotal;

            if Ni < phaseFrac(1)
                skipMask(ki) = true;
                continue;
            end

            [uniquePh, ~, ic] = unique(validPh);
            cnts = accumarray(ic, 1);
            [maxCnt, idxMax] = max(cnts);
            domPhase = uniquePh(idxMax);
            fracDom = maxCnt / numel(validPh);

            if fracDom < phaseFrac(2)
                skipMask(ki) = true;
                continue;
            end

            oriList = reshape(winOri, [], 4);

            validFull = kernelFull & winSample & ~winProt;
            validFullVec = validFull(:);

            phaseVec = winPhase(:);
            oriSub = oriList(validFullVec, :);
            phaseSub = phaseVec(validFullVec);

            domMask = (phaseSub == domPhase);
            nqList = oriSub(domMask, :);

            if size(nqList, 1) < 2
                skipMask(ki) = true;
                continue;
            end

            qMean = calc_mean_ori_hole_strict( ...
                nqList, winPhase, domPhase, misTol, innerK, validOuter, r, validFullVec, threshFrac);

            if isempty(qMean)
                skipMask(ki) = true;
                continue;
            end

            localPhase(ki) = domPhase;
            localOri(ki, :) = [qMean.a qMean.b qMean.c qMean.d];
            localFilled(ki) = true;

            patchPhase(iL, jL) = domPhase;
            patchOri(iL, jL, :) = reshape([qMean.a qMean.b qMean.c qMean.d], 1, 1, 4);

            madeChange = true;

            for kj = 1:nPix
                if ~skipMask(kj)
                    continue;
                end
                if abs(iL_vec(kj) - iL) <= 2*r && abs(jL_vec(kj) - jL) <= 2*r
                    skipMask(kj) = false;
                end
            end
        end
    end

    result.pixList = pixList;
    result.localPhase = localPhase;
    result.localOri = localOri;
    result.localFilled = localFilled;
end

function [meanOri, clusterSizes] = calc_mean_ori_hole_strict(nqList, winPhase, domPhase, misTol, innerK, outerMaskForRing, r, validFullVec, threshFrac)
    % Purpose:
    %   Computes a local representative orientation for hole filling.
    %
    % Logic:
    %   Dominant-phase neighbour orientations are clustered. If one cluster is
    %   sufficiently dominant it is used directly; otherwise a ring-based
    %   fallback attempts to recover a coherent local orientation.
    Ndom = size(nqList, 1);

    if Ndom == 0
        meanOri = [];
        clusterSizes = [];
        return;
    end

    if Ndom == 1
        meanOri = quaternion(nqList(1,1), nqList(1,2), nqList(1,3), nqList(1,4));
        clusterSizes = 1;
        return;
    end

    D = 2 * acos(min(abs(nqList * nqList.'), 1));
    D(1:size(D,1)+1:end) = 0;

    Z = linkage(squareform(D), 'single');
    idx = cluster(Z, 'cutoff', misTol, 'criterion', 'distance');

    counts = accumarray(idx, 1);
    [maxCount, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;

    if (maxCount / Ndom) >= threshFrac
        meanOri = mean(quaternion(nqList(members, :).'), 'meanOrientation');
        return;
    end

    distMap = round(bwdist(innerK));
    distMap(~outerMaskForRing) = 0;

    validLinear = find(validFullVec & (winPhase(:) == domPhase));

    for ringW = 1:(r-1)
        ringMask = (distMap == ringW);
        ringLinear = find(ringMask(:));
        [~, commonIdx] = intersect(validLinear, ringLinear);

        if numel(commonIdx) < 2
            continue;
        end

        if max(commonIdx) > size(nqList, 1)
            continue;
        end

        qRing = nqList(commonIdx, :);

        Dr = 2 * acos(min(abs(qRing * qRing.'), 1));
        Dr(1:size(Dr,1)+1:end) = 0;

        Zr = linkage(squareform(Dr), 'single');
        ir = cluster(Zr, 'cutoff', misTol, 'criterion', 'distance');
        cr = accumarray(ir, 1);
        [mxR, dcR] = max(cr);

        if (mxR / numel(commonIdx)) >= 0.5
            mbR = (ir == dcR);
            meanOri = mean(quaternion(qRing(mbR, :).'), 'meanOrientation');
            return;
        end
    end

    meanOri = [];
end
%% =========================================================================
%  HOLE FILLING - RELAXED (MPF)
%% =========================================================================
function [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFillingMPF(ebsd, oriQuatClean, phaseMapClean, protectedMask, sampleMask)
    % Purpose:
    %   Fills unprotected notIndexed holes using a serial MPF strategy with
    %   vectorised candidate screening and local orientation validation.
    %
    % Inputs:
    %   ebsd          - EBSD map
    %   oriQuatClean  - quaternion orientation map
    %   phaseMapClean - cleaned phase map
    %   protectedMask - pixels excluded from ordinary filling
    %   sampleMask    - valid specimen mask
    %
    % Outputs:
    %   ebsd_fill     - EBSD map after BFS filling
    %   phaseMapClean - updated phase map
    %   oriQuatClean  - updated quaternion map
    %
    % Logic:
    %   Candidate holes are screened using vectorised neighbourhood support
    %   maps. Only screened survivors undergo local dominant-phase checking and
    %   orientation-mean estimation. The padded maps are updated in-pass so that
    %   newly filled pixels can support later decisions within the same pass.
    global params

    t0 = tic;
    fprintf('\nHole Filling MPF (serial, vectorised screening):\n');

    scaleLead = params.scaleLead;
    minLead   = params.minLead;
    misTol_ori = params.misTol_ori;

    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);

    [Nrow, Ncol] = size(phaseMapClean);
    radii = params.radius_fill;

    ebsdPhaseIdGrid = reshape(ebsd.phaseId, Nrow, Ncol);
    initialHoleMask = (ebsdPhaseIdGrid == notIndexedId) & ~protectedMask & sampleMask;

    for radius = radii
        tRad = tic;
        fprintf('--- MPF (radius = %d) ---\n', radius);

        phaseFrac = params.phaseFrac(radius);
        N = 2*radius + 1;

        kernelFull = double(fspecial('disk', radius)) > 0;
        innerKernel = false(N);
        innerKernel(radius+1, radius+1) = true;
        outerKernel = kernelFull & ~innerKernel;
        nOuter = sum(outerKernel(:));

        % Build padded maps ONCE per radius, then update in place
        paddedPhase = padarray(phaseMapClean, [radius radius], 0, 'both');
        paddedOri   = padarray(oriQuatClean, [radius radius], 0, 'both');

        passes = 0;
        filledLastPass = false(Nrow, Ncol);
        fillCountRadius = 0;

        while true
            passes = passes + 1;
            fprintf('\n--- Hole Filling Pass %d (Radius = %d) ---\n', passes, radius);

            baseHoleMask = (phaseMapClean == notIndexedId) & ~protectedMask & sampleMask;

            if ~any(baseHoleMask(:))
                fprintf('⚠ No more holes remain at radius = %d\n', radius);
                fprintf('✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
                break;
            end

            if passes == 1
                holeMask = baseHoleMask;
            else
                dilatedFilled = imdilate(filledLastPass, outerKernel);
                holeMask = baseHoleMask & dilatedFilled;
            end

            if ~any(holeMask(:))
                fprintf('⚠ No more discoverable holes at radius = %d\n', radius);
                fprintf('✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
                break;
            end

            % --------------------------------------------------------------
            % Vectorised cheap screening for this pass
            % --------------------------------------------------------------
            validIndexedMask = (phaseMapClean > 0) & (phaseMapClean ~= notIndexedId) & sampleMask & ~protectedMask;
            neighCount = conv2(double(validIndexedMask), double(outerKernel), 'same');
            Ni_map = neighCount / nOuter;

            numPhases = numel(MinPhaseIds);
            phaseCountMaps = zeros(Nrow, Ncol, numPhases, 'double');

            for pp = 1:numPhases
                pid = MinPhaseIds(pp);
                phaseMask = (phaseMapClean == pid) & sampleMask & ~protectedMask;
                phaseCountMaps(:,:,pp) = conv2(double(phaseMask), double(outerKernel), 'same');
            end

            [maxPhaseCount, domIdx] = max(phaseCountMaps, [], 3);
            domPhaseMap = zeros(Nrow, Ncol);
            for pp = 1:numPhases
                domPhaseMap(domIdx == pp) = MinPhaseIds(pp);
            end
            fracDomMap = maxPhaseCount ./ max(neighCount, 1);

            candidateMask = holeMask & ...
                            (Ni_map >= phaseFrac(1)) & ...
                            (fracDomMap >= phaseFrac(2));

            [candRows, candCols] = find(candidateMask);
            numCandidates = numel(candRows);

            if numCandidates == 0
                fprintf('⚠ No screened candidates at radius = %d, pass = %d\n', radius, passes);
                fprintf('✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
                break;
            end

            fprintf('Total screened candidates: %d\n', numCandidates);

            fillCountPass = 0;
            filledLastPass = false(Nrow, Ncol);
            tLoop = tic;
            nextPrint = 1;
            chunkStart = 1;
            chunkFilled = 0;

            % --------------------------------------------------------------
            % Serial fill only on screened survivors
            % --------------------------------------------------------------
            for k = 1:numCandidates
                i = candRows(k);
                j = candCols(k);

                % Pixel may already have been filled earlier this pass
                if phaseMapClean(i,j) ~= notIndexedId
                    continue;
                end

                iP = i + radius;
                jP = j + radius;

                winPhase = paddedPhase(iP-radius:iP+radius, jP-radius:jP+radius);
                winOri   = paddedOri(iP-radius:iP+radius, jP-radius:jP+radius, :);

                maskValid = (winPhase > 0 & winPhase ~= notIndexedId) & outerKernel;
                neighPhases = winPhase(maskValid);

                if isempty(neighPhases)
                    continue;
                end

                Ni = numel(neighPhases) / nOuter;
                if Ni < phaseFrac(1)
                    continue;
                end

                % Dominant phase from current local window
                % (recomputed locally so in-pass updates can influence outcome)
                localCounts = zeros(numPhases, 1);
                for pp = 1:numPhases
                    localCounts(pp) = sum(neighPhases == MinPhaseIds(pp));
                end
                [maxCount, idxMax] = max(localCounts);

                if maxCount == 0
                    continue;
                end

                domPhase = MinPhaseIds(idxMax);
                fracDom = maxCount / numel(neighPhases);

                if fracDom < phaseFrac(2)
                    continue;
                end

                % Early cheap skip before expensive mean calculation
                if maxCount < 2
                    continue;
                end

                oriList = reshape(winOri, [], 4);
                qList = oriList(maskValid(:), :);

                qMean = calc_mean_ori_hole_relaxed( ...
                    qList, neighPhases, domPhase, misTol_ori, ...
                    innerKernel, outerKernel, radius, maskValid, ...
                    scaleLead, minLead);

                if isempty(qMean)
                    pct = floor(100 * k / numCandidates);
                    if pct >= nextPrint
                        elapsed = toc(tLoop);
                        scannedChunk = k - chunkStart + 1;
                        fprintf('      %3d%% | scanned = %d | filled = %d | chunk fill = %.1f%% | elapsed = %.1f s\n', ...
                            pct, scannedChunk, chunkFilled, ...
                            100 * chunkFilled / max(1, scannedChunk), elapsed);
                        nextPrint = pct + 1;
                        chunkStart = k + 1;
                        chunkFilled = 0;
                    end
                    continue;
                end

                phaseMapClean(i,j) = domPhase;
                oriQuatClean(i,j,:) = [qMean.a qMean.b qMean.c qMean.d];

                % Update padded maps immediately so later checks in this pass
                % can benefit from newly filled pixels
                paddedPhase(iP,jP) = domPhase;
                paddedOri(iP,jP,:) = [qMean.a qMean.b qMean.c qMean.d];

                filledLastPass(i,j) = true;
                fillCountPass = fillCountPass + 1;
                fillCountRadius = fillCountRadius + 1;
                chunkFilled = chunkFilled + 1;

                pct = floor(100 * k / numCandidates);
                if pct >= nextPrint
                    elapsed = toc(tLoop);
                    scannedChunk = k - chunkStart + 1;
                    fprintf('      %3d%% | scanned = %d | filled = %d | chunk fill = %.1f%% | elapsed = %.1f s\n', ...
                        pct, scannedChunk, chunkFilled, ...
                        100 * chunkFilled / max(1, scannedChunk), elapsed);
                    nextPrint = pct + 1;
                    chunkStart = k + 1;
                    chunkFilled = 0;
                end
            end

            fprintf('✔ Pass %d complete: filled %d holes\n', passes, fillCountPass);

            if fillCountPass == 0
                fprintf('⚠ None of the discovered holes filled.\n');
                fprintf('✔ Total holes filled by radius %d: %d\n', radius, fillCountRadius);
                break;
            end
        end

        fprintf('✔ Radius %d total complete: filled %d holes (%.2f s)\n', ...
            radius, fillCountRadius, toc(tRad));
    end

    % ---------------------------------------------------------------------
    % Rebuild EBSD object
    % ---------------------------------------------------------------------
    ebsd_fill = ebsd;
    ebsd_fill.phaseId(:) = phaseMapClean(:);

    qFull_fill = quaternion( ...
        oriQuatClean(:,:,1), ...
        oriQuatClean(:,:,2), ...
        oriQuatClean(:,:,3), ...
        oriQuatClean(:,:,4));

    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (ebsd_fill.phaseId == pid) & initialHoleMask(:);
        if ~any(mask)
            continue;
        end
        ebsd_fill(mask).orientations = orientation(qFull_fill(mask), ebsd(pname).CS);
    end

    fprintf('✔ MPF complete (%.2f s total)\n', toc(t0));
end




function [meanOri, clusterSizes] = calc_mean_ori_hole_relaxed(qList, neighPhases, domPhase, misTol_ori, innerKernel, outerKernel, radius, maskValid, scaleLead, minLead)
    % Purpose:
    %   Computes a local mean orientation for relaxed MPF hole filling.
    %
    % Logic:
    %   Dominant-phase neighbours are clustered first globally, then optionally
    %   re-evaluated using ring-wise support if no cluster has a sufficient lead.
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
    D = 2 * acos(min(abs(qListDom * qListDom.'), 1));
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);
    Z = linkage(Dc, 'single');
    idx = cluster(Z, 'cutoff', misTol_ori, 'criterion', 'distance');

    % --- Step 3: Dominant cluster check (Lead) ---
    counts = accumarray(idx, 1);
    [maxCount, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;

    otherCounts = counts;
    otherCounts(domCluster) = 0;
    leadDiff = maxCount - max(otherCounts, 0);

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
        if ~any(ringMask(:))
            continue;
        end

        validLinearIdx = find(maskValid);
        ringLinearIdx = find(ringMask(:));
        [~, commonidx, ~] = intersect(validLinearIdx, ringLinearIdx);

        qRing = qList(commonidx, :);
        neighPhasesRing = neighPhases(commonidx);
        qRingDom = qRing(neighPhasesRing == domPhase, :);
        Nring = size(qRingDom, 1);

        if Nring < 2
            continue;
        end

        D_ring = 2 * acos(min(abs(qRingDom * qRingDom.'), 1));
        D_ring(1:size(D_ring,1)+1:end) = 0;
        Dc_ring = squareform(D_ring);
        Z_ring = linkage(Dc_ring, 'single');
        idx_ring = cluster(Z_ring, 'cutoff', misTol_ori, 'criterion', 'distance');
        counts_ring = accumarray(idx_ring, 1);
        [maxCountRing, domClusterRing] = max(counts_ring);

        other_ring = counts_ring;
        other_ring(domClusterRing) = 0;
        leadDiffRing = maxCountRing - max(other_ring, 0);
        reqLeadRing = max(minLead, ceil(scaleLead * Nring));

        if leadDiffRing >= reqLeadRing
            membersRing = (idx_ring == domClusterRing);
            qClusterRing = quaternion(qRingDom(membersRing,:).');
            meanOri = mean(qClusterRing, 'meanOrientation');
            found = true;
            break;
        end
    end

    % --- Step 5: fallback ---
    if ~found
        meanOri = [];
    end
end

function showPhaseStats(ebsdObj, phases, msg)
    % Purpose:
    %   Prints phase fractions and counts for the supplied EBSD object.
    fprintf('\n%s\n', msg);
    fprintf('--------------------------------\n');
    total = numel(ebsdObj);
    for i = 1:numel(phases)
        n = numel(ebsdObj(phases{i}));
        fprintf('%-12s: %6d points (%.2f%%)\n', phases{i}, n, 100*n/total);
    end
    fprintf('--------------------------------\n');
end

function plotPhaseMap(ebsdObj, sampleName, exportPath, suffix, res)
    % Purpose:
    %   Exports a phase map figure for the current cleaning stage.
    f = figure('Visible','off');
    plot(ebsdObj,'phase');
    leg = legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','on','FontSize',10);
    leg.Position(1) = 0.5 - leg.Position(3)/2;
    savePNG(f, sprintf('%s_PhaseMap_%s',sampleName,suffix), exportPath, res);
end

function plotIPFMapPhases(ebsdObj, sampleName, exportPath, suffix, res)
    % Purpose:
    %   Exports per-phase IPF maps for the current cleaning stage.
    phases = ebsdObj.mineralList;
    for i = 1:numel(phases)
        pname = phases{i};
        if strcmpi(pname,'notIndexed'), continue; end
        ebsdPhase = ebsdObj(pname);
        if isempty(ebsdPhase), continue; end
        f = figure('Visible','off');
        plot(ebsdPhase, ebsdPhase.orientations);
        axis equal;
        savePNG(f, sprintf('%s_IPFMap_%s_%s', sampleName, pname, suffix), exportPath, res);
    end
end

function savePNG(figHandle, filenameStem, exportPath, res)
    % Purpose:
    %   Saves a figure as a PNG at the requested export resolution.
    exportgraphics(figHandle, fullfile(exportPath,[filenameStem '.png']), 'Resolution', res);
    close(figHandle);
    fprintf('Saved: %s.png\n', filenameStem);
end
