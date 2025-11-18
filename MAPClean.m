%% Modular EBSD Cleaning Pipeline
clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

dataDir = fullfile(pwd,'DataFiles');
fileList = dir(fullfile(dataDir,'01a6.ctf'));
checkpointDir = fullfile(pwd,'checkpoints');
exportDir     = fullfile(pwd,'exports');
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end
if ~exist(exportDir,'dir'), mkdir(exportDir); end

%% --- Stage control ---
runStart    = false;
runMAD      = false;
runPhaseWSR = false;
runOriWSR   = false;
runHoleFill = true;
runSaveFile = true;

%% --- Parameters ---
global params
params.thresholdFrac = 0.75; % proportion of single phase to be considered dominant for mean orientation
params.exportRes    = 300;   % export resolution in dpi
params.madThreshold = 0.9;   % radians
params.radius_phase = 2;     % radius of grid to detect phase WSR
params.radius_ori   = 2;     % radius of grid to detect orientation WSR
params.misTol_ori   = 5*degree;
params.minFrac_ori  = 0.25;
params.radius_fill = [6 5 4 3 2 1]; % Maximum neighbourhood radius for hole filling
params.min_neighbours = 3;     
params.min_dom_frac   = 0.50;  

% Adaptive phase fraction based on radius
params.phaseFrac = containers.Map('KeyType','double','ValueType','any');
params.phaseFrac(6) = [0.4 0.75];
params.phaseFrac(5) = [0.4 0.75];
params.phaseFrac(4) = [0.4 0.75];
params.phaseFrac(3) = [0.4 0.75];
params.phaseFrac(2) = [0.4 0.75];
params.phaseFrac(1) = [0.4 0.75];

disp('✔ Parameters initialised');

%% --- Loop over samples ---
for fi = 1:numel(fileList)
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    exportPath = fullfile(exportDir, sampleName);
    if ~exist("exportPath",'dir'), mkdir(exportPath); end
    diaryFile = fullfile(exportDir, [sampleName '_logfile.txt']);
    diary(diaryFile)
    diary on
    fprintf('\n===== Sample: %s =====\n', sampleName);

    ebsd_raw = EBSD.load(fullfile(dataDir,fileList(fi).name),'convertSpatial2EulerReferenceFrame').gridify;
    phases = ebsd_raw.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);

    % Set Phase Colours
    ForsteriteColor = str2rgb('red');
    DiopsideColor   = str2rgb('blue');
    AnorthiteColor  = str2rgb('yellow');
    for i = 1:numel(MinPhaseNames)
        pname = MinPhaseNames{i};
        if strcmpi(pname, 'Forsterite')
            ebsd_raw(pname).color = ForsteriteColor;
        elseif strcmpi(pname, 'Diopside')
            ebsd_raw(pname).color = DiopsideColor;
        elseif strcmpi(pname, 'Anorthite')
            ebsd_raw(pname).color = AnorthiteColor;
        end
    end
    [Nrow, Ncol] = ebsd_raw.size;
    clear AnorthiteColor DiopsideColor ForsteriteColor

    % --- Checkpoint file paths ---
    madFile   = fullfile(checkpointDir, sprintf('%s_ebsd_mad.mat', sampleName));
    phaseFile = fullfile(checkpointDir, sprintf('%s_ebsd_phase.mat', sampleName));
    oriFile   = fullfile(checkpointDir, sprintf('%s_ebsd_ori.mat', sampleName));
    fillFile  = fullfile(checkpointDir, sprintf('%s_ebsd_fill.mat', sampleName));
    paramFile = fullfile(checkpointDir, sprintf('%s_params.mat', sampleName));

    ebsd = ebsd_raw;

    %% --- Starting summary and plots ---
    if runStart
        showPhaseStats(ebsd_raw, phases, 'Phase distribution before cleaning');
        plotPhaseMap(ebsd_raw, sampleName, exportPath, 'raw', params.exportRes);
        plotIPFMapPhases(ebsd_raw, sampleName, exportPath, 'raw', params.exportRes);
    end

    %% --- MAD Filter ---
    if runMAD
        [ebsd_mad, badPixels] = doMADFilter(ebsd, sampleName, exportPath);
        save(madFile,'ebsd_mad');
        fprintf('✔ Checkpoint MAD saved\n');
        ebsd = ebsd_mad;
    elseif exist(madFile,'file')
        load(madFile,'ebsd_mad');
        ebsd = ebsd_mad;
        fprintf('\nMAD checkpoint loaded.\n');
    else
        disp('MAD filter skipped and no previous file exists.');
    end

    %% --- Check data quality ---
    phases = ebsd.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    fracNotIndexed = sum(ebsd.phaseId == notIndexedId) / numel(ebsd.phaseId);

    fprintf('\n--- Data Quality Assessment ---\n');
    if fracNotIndexed < 0.75
        runStrict = true;
        fprintf('✔ Data sufficiently indexed (%.2f%% notIndexed). STRICT mode selected.\n', fracNotIndexed*100);
    else
        runStrict = false;
        fprintf('⚠ Sparse data (%.2f%% notIndexed). RELAXED mode selected.\n', fracNotIndexed*100);
    end

    %% --- WSR Processing ---
    fprintf('\n=== WSR Processing ===\n');
    if runStrict
        % --- Phase WSR ---
        if runPhaseWSR
            fprintf('\n=== Running Phase WSR ===\n');
            [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR(ebsd, sampleName, exportPath);
            save(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            fprintf('✔ Phase WSR saved\n');
            ebsd = ebsd_phase;
        elseif exist(phaseFile,'file')
            load(phaseFile,'ebsd_phase','phaseMapClean','oriQuatClean','protectedMask');
            fprintf('✔ Loaded existing Phase WSR checkpoint.\n');
            ebsd = ebsd_phase;
        else
            fprintf('⚠ Phase WSR skipped. Using MAD output.\n');
        end

        % --- Orientation WSR ---
        if runOriWSR
            if ~exist('ebsd_phase','var') && ~exist(phaseFile,'file')
                fprintf('⚠ Orientation WSR cannot run without Phase WSR. Skipping.\n');
            else
                fprintf('\n=== Running Orientation WSR ===\n');
                [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath);
                save(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean');
                fprintf('✔ Orientation WSR saved\n');
                ebsd = ebsd_ori;
            end
        elseif exist(oriFile,'file')
            load(oriFile,'ebsd_ori','phaseMapClean','oriQuatClean');
            fprintf('✔ Loaded existing Orientation WSR checkpoint.\n');
            ebsd = ebsd_ori;
        else
            fprintf('⚠ Orientation WSR skipped.\n');
        end

    else
        fprintf('⚠ No WSR in RELAXED mode. Skipping WSR entirely.\n');
        [phaseMapClean, oriQuatClean, protectedMask] = residuals(ebsd);
    end

    %% --- Hole Filling ---
    fprintf('\n=== Filling Holes ===\n');
    if runHoleFill
        ebsd_fill = ebsd;
        [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFilling(ebsd_fill, oriQuatClean, phaseMapClean, params.radius_fill, protectedMask);
        ebsd_fill = EBSD(ebsd_fill, 'convert');
        showPhaseStats(ebsd_fill, phases, 'Phase distribution after hole fill');
        plotPhaseMap(ebsd_fill, sampleName, exportPath, 'Fill', params.exportRes);
        plotIPFMapPhases(ebsd_fill, sampleName, exportPath, 'Fill', params.exportRes);
        save(fillFile,'ebsd_fill','phaseMapClean','oriQuatClean','protectedMask');
        fprintf('✔ Checkpoint Hole Fill saved\n');
        ebsd = ebsd_fill;

    elseif exist(fillFile,'file')
        load(fillFile,'ebsd_fill','phaseMapClean','oriQuatClean','protectedMask');
        ebsd = ebsd_fill;
        fprintf('\nHole Fill checkpoint loaded.\n');
    else
        disp('Hole Filling skipped and no previous file exists.');
    end

    %% --- Export final cleaned EBSD ---
    if runSaveFile
        outFile = fullfile(dataDir,[sampleName '_clean.ctf']);
        ebsd.export(outFile);
        fprintf('✔ Saved cleaned EBSD: %s\n', outFile);
        save(paramFile,'params');
        fprintf('✔ Parameters saved in a mat file\n');
    end
    diary off
end


%% MAIN FUNCTIONS
% ---------- MAD filter ----------
function [ebsd_mad, badPixels] = doMADFilter(ebsd, sampleName, exportPath)
    global params;
    fprintf('\n=== Applying MAD Filter (Threshold = %.2f rad) ===\n', params.madThreshold);
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

% ---------- Phase Wild Spike Removal ----------
function [ebsd_phase, phaseMapClean, oriQuatClean, protectedMask] = doPhaseWSR(ebsd, sampleName, exportPath)
    global params;
    fprintf('\n=== Starting Phase WSR (Radius = %d) ===\n', params.radius_phase);
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
        mask = (phaseMapoG == pid);
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
        mask = (phaseMapoG == pid);
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
                        [meanOri, ~] = calc_mean_ori_wsr(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
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
                    [meanOri, ~] = calc_mean_ori_wsr(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
                    oriQuatClean(i,j,:) = [meanOri.a, meanOri.b, meanOri.c, meanOri.d];
                end
            end
        end
    end
    % Update EBSD object
    ebsd_phase = ebsd;
    ebsd_phase.phaseId(:) = phaseMapClean(:);
    qFull_phase = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        mask = (phaseMapClean==pid);
        if ~any(mask,'all'), continue; end
        ebsd_phase(mask).orientations = orientation(qFull_phase(mask), ebsd(pname).CS);
    end
    % Show stats and plots
    showPhaseStats(ebsd_phase, phases, 'Phase distribution after Phase WSR');
    plotPhaseMap(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);
    plotIPFMapPhases(ebsd_phase, sampleName, exportPath, 'PhaseWSR', params.exportRes);   
end

% ---------- Orientation Wild Spike Removal ----------
function [ebsd_ori, oriQuatClean, phaseMapClean] = doOrientationWSR(ebsd, oriQuatClean, phaseMapClean, MinPhaseIds, sampleName, exportPath)  
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
        mask = (phaseMapClean == pid);
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
                end
                continue;
            end
            % --- N >= 3 neighbours ---
            dots = abs(neighbourQuatsList * currentQ_vec'); dots(dots>1)=1;
            misAngles = 2*acos(dots);
            if all(misAngles < misTol_ori)
                qMean = mean(quaternion(neighbourQuatsList.'), 'meanOrientation');
            else
                qMean = calc_mean_ori_wsr(neighbourQuatsList, params.misTol_ori, Nneighbours, currentQ_vec);
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
                if isTwin continue; end % Skip if it’s a twin
            end
            % --- Update orientation ---
            oriQuatClean(i,j,:) = [qMean.a qMean.b qMean.c qMean.d];
            paddedOri(iP,jP,:)  = [qMean.a qMean.b qMean.c qMean.d];
        end
    end
    % --- Push back to EBSD ---
    ebsd_ori = ebsd;
    ebsd_ori.phaseId(:) = phaseMapClean(:);
    qFull_wsr = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p); pname = MinPhaseNames{p};
        mask = (phaseMapClean==pid);
        if ~any(mask,'all'), continue; end
        ebsd_ori(mask).orientations = orientation(qFull_wsr(mask), ebsd_ori(pname).CS);
    end
    % Show stats and plots
    showPhaseStats(ebsd_ori, phases, 'Phase distribution after Ori WSR');
    plotPhaseMap(ebsd_ori, sampleName, exportPath, 'OriWSR', params.exportRes);
    plotIPFMapPhases(ebsd_ori, sampleName, exportPath, 'OriWSR', params.exportRes);
end

%% ---------- Hole Filling ----------
function [ebsd_fill, phaseMapClean, oriQuatClean] = doHoleFilling(ebsd, oriQuatClean, phaseMapClean, radii, protectedMask)
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
        % initial hole mask: holes that are eligible for processing (not protected)
        holeMask = (phaseMapClean == notIndexedId) & ~protectedMask;
        clusterId = 0;
        % While there exists an undiscovered hole, form a cluster and fill it
        while true
            % find first undiscovered hole
            [r0, c0] = find(holeMask & ~visited, 1, 'first');
            if isempty(r0)
                break; % no more undiscovered holes at this radius
                fprintf('\n No more undiscovered holes in radius = %d \n', radius);
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
                        if ni < 1 || nj < 1 || ni > Nrow || nj > Ncol, continue; end
                        % neighbour is a hole (notIndexed) and not protected
                        if holeMask(ni,nj) && ~visited(ni,nj)
                            visited(ni,nj) = true;
                            queue(end+1, :) = [ni, nj];
                        end
                    end
                end
            end
            % queue now contains the full connected cluster (all hole pixels connected to start)
            totalCandidates = size(queue,1);

            % ---------- FILL PHASE (attempt to fill each pixel in the cluster) ----------
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
                % Align kernel to window: compute kernel slice that overlaps window
                kRowStart = 1 + (rmin - (i - radius)); % when rmin == i-radius -> 1
                kColStart = 1 + (cmin - (j - radius));
                kRowEnd   = kRowStart + (rmax - rmin);
                kColEnd   = kColStart + (cmax - cmin);
                kernelLocal = outerKernel(kRowStart:kRowEnd, kColStart:kColEnd);
                % Valid neighbor mask inside window: indexed and not 'notIndexed'
                validMask = (winPhase > 0 & winPhase ~= notIndexedId) & kernelLocal;
                neighPhases = winPhase(validMask);

                % Condition 1: check fraction of valid neighbours (Ni threshold)
                Ni = numel(neighPhases) / max(1, numNeighbours); % avoid divide by zero
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

                % Passed thresholds → assign phase
                phaseMapClean(i,j) = domPhase;
                filledCount = filledCount + 1;

                % Compute mean orientation using q's of valid dominated-phase neighbours
                qList = reshape(winOri, [], 4); % rows correspond to same ordering as winPhase(:)
                currentQ_vec = squeeze(oriQuatClean(i,j,:))';

                % Use existing helper (calc_mean_ori_hole) to compute qMean (safe if qValid empty / small)
                qMean = calc_mean_ori_hole(qList, neighPhases, domPhase, misTol_ori, currentQ_vec, innerKernel, outerKernel, radius, validMask);
                oriQuatClean(i,j,:) = [qMean.a, qMean.b, qMean.c, qMean.d];
            end

            % Print cluster result if filled anything
            if filledCount > 0 & totalCandidates > 10
                fprintf('Cluster %d: filled %d/%d\n', clusterId, filledCount, totalCandidates);
            end
            % Update holeMask for next cluster (we exclude protected pixels and those already filled)
            holeMask = (phaseMapClean == notIndexedId) & ~protectedMask;
            % Note: visited remains true for discovered cluster pixels so we won't rediscover same cluster
        end 
    end
    % Rebuild EBSD object with final maps
    ebsd_fill = ebsd;
    ebsd_fill.phaseId(:) = phaseMapClean(:);
    qFull_fill = quaternion(oriQuatClean(:,:,1), oriQuatClean(:,:,2), oriQuatClean(:,:,3), oriQuatClean(:,:,4));
    for p = 1:numel(MinPhaseIds)
        pid = MinPhaseIds(p);
        pname = MinPhaseNames{p};
        mask = (phaseMapClean == pid);
        if ~any(mask,'all'), continue; end
        ebsd_fill(mask).orientations = orientation(qFull_fill(mask), ebsd(pname).CS);
    end
end


%% =================== Helper Functions ===================
function showPhaseStats(ebsdObj, phases, msg)
    fprintf('\n%s\n', msg); fprintf('--------------------------------\n');
    total = numel(ebsdObj);
    for i=1:numel(phases)
        n=numel(ebsdObj(phases{i}));
        fprintf('%-12s: %6d points (%.2f%%)\n', phases{i}, n, 100*n/total);
    end
    fprintf('--------------------------------\n');
end
function plotPhaseMap(ebsdObj, sampleName, exportPath, suffix,res)
    f=figure('Visible','off'); plot(ebsdObj,'phase');
    leg=legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','on','FontSize',10);
    leg.Position(1)=0.5-leg.Position(3)/2;
    savePNG(f,sprintf('%s_PhaseMap_%s',sampleName,suffix),exportPath,res);
end
function plotIPFMapPhases(ebsdObj, sampleName, exportPath, suffix, res)
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
    exportgraphics(figHandle,fullfile(exportPath,[filenameStem,'.png']),'Resolution',res);
    close(figHandle);
    fprintf('Saved: %s.png\n', filenameStem);
end

function [meanOri, clusterSizes] = calc_mean_ori_wsr(qList, misTol_ori, Nneighbours, currentQ_vec)
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
    [~, domCluster] = max(counts);
    members = (idx == domCluster);
    clusterSizes = counts;
    % Step 5: compute representative orientation based on dominant cluster fraction
    domClusterFrac = counts(domCluster) / Nneighbours;
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
%%
function [meanOri, clusterSizes] = calc_mean_ori_hole(qList, neighPhases, domPhase, misTol_ori, currentQ_vec, innerKernel, outerKernel, radius, maskValid)
    global params
    thresholdFrac = params.thresholdFrac;

    % 1. Filter to only dominant phase neighbours
    qListDom = qList(neighPhases == domPhase, :);
    Ndom = size(qListDom,1);

    % Convert current orientation to quaternion
    currentQ = quaternion(currentQ_vec(1), currentQ_vec(2), currentQ_vec(3), currentQ_vec(4));

    if Ndom == 0
        meanOri = currentQ;
        clusterSizes = 0;
        return;
    end

    % --- Step 1: pairwise angular distances ---
    D = 2 * acos(min(abs(qListDom*qListDom.'),1));
    D(1:size(D,1)+1:end) = 0;
    Dc = squareform(D);

    % --- Step 2: hierarchical clustering ---
    Z = linkage(Dc,'single');

    % --- Step 3: cluster neighbours by misorientation tolerance ---
    idx = cluster(Z,'cutoff',misTol_ori,'criterion','distance');

    % --- Step 4: select dominant cluster ---
    counts = accumarray(idx,1);
    [maxCount, domCluster] = max(counts);
    members = (idx==domCluster);
    clusterSizes = counts;

    % Step 5: compute representative orientation
    domFrac = maxCount/Ndom;

    if domFrac >= thresholdFrac
        qCluster = quaternion(qListDom(members,:).'); 
        meanOri = mean(qCluster, 'meanOrientation');
        return;
    end

    % --- Ring-based fallback if dominant cluster too small ---
    distMap = bwdist(innerKernel);
    distMap = round(distMap);
    distMap(~outerKernel) = 0;

    found = false;
    for ringWidth = 2:(radius-1)
        ringMask = (distMap==ringWidth);
        if ~any(ringMask(:)), continue; end
        ringMaskVec = ringMask(:);
        validLinearIdx = find(maskValid);
        ringLinearIdx = find(ringMaskVec);
        % Find indices that are both valid and in the ring
        [~, commonidx] = intersect(validLinearIdx, ringLinearIdx); % indices in validLinearIdx that are also in ringLinearIdx.
        % Now extract qList and neighPhases safely
        qRing = qList(commonidx, :);
        neighPhasesRing = neighPhases(commonidx);
        % Keep only dominant phase in this ring
        qRing = qRing(neighPhasesRing==domPhase,:);
        Nring = size(qRing,1);
        if Nring<2, continue; end

        D_ring = 2*acos(min(abs(qRing*qRing.'),1));
        D_ring(1:size(D_ring,1)+1:end) = 0;
        Dc_ring = squareform(D_ring);
        Z_ring = linkage(Dc_ring,'single');
        idx_ring = cluster(Z_ring,'cutoff',misTol_ori,'criterion','distance');
        counts_ring = accumarray(idx_ring,1);
        [maxCountRing, domClusterRing] = max(counts_ring);
        domFracRing = maxCountRing/Nring;
        if domFracRing>=0.5
            membersRing = (idx_ring==domClusterRing);
            qClusterRing = quaternion(qRing(membersRing,:).');
            meanOri = mean(qClusterRing,'meanOrientation');
            found = true;
            break;
        end
    end

    if ~found
        meanOri = currentQ;
    end
end

%%
function [phaseMapClean, oriQuatClean, protectedMask] = residuals(ebsd)
    global params
    fprintf('Initialising phaseMapClean and oriQuatClean from EBSD.\n');
    phases = ebsd.mineralList;
    [Nrow, Ncol] = ebsd.size;
    phaseMapClean = reshape(double(ebsd.phaseId), Nrow, Ncol);
    oriQuatClean = zeros(Nrow, Ncol, 4);
    N_total = Nrow * Ncol;
    for p = 1:numel(phases)
        pname = phases{p};
        mask = ebsd(pname).isIndexed;
        if ~any(mask), continue; end
        q = quaternion(ebsd(pname).orientations);
        idx = find(mask);
        oriQuatClean(idx)             = [q.a];
        oriQuatClean(idx + N_total)   = [q.b];
        oriQuatClean(idx + 2*N_total) = [q.c];
        oriQuatClean(idx + 3*N_total) = [q.d];
    end
    protectedMask = ebsd.mad > params.madThreshold;  
end

%%
