%% GasLeakBot


clear; clc; close all;
rng(1);                              % fixed random seed -> same noise every run (repeatable)

%% 1) Parameters

% ---- World and grid ----
res       = 0.25;                    % size of one grid cell (m)
worldSize = 20;                      % world is 20 m x 20 m
N         = worldSize/res;           % cells per side (80)
dt        = 0.1;                     % simulation time step (s)
maxSteps  = 6000;                    % safety limit

% ---- LiDAR ----
nRays  = 72;                         % number of beams
maxR   = 8;                          % maximum range (m)
relAng = linspace(-pi, pi, nRays+1); % beam angles relative to the robot heading
relAng(end) = [];                    % remove duplicate (-pi and +pi are the same beam)

% ---- Robot motion and obstacle avoidance ----
start   = [2.5; 2; pi/4];            % [x; y; heading]  (heading in radians)
vMax    = 0.2;                       % maximum forward speed (m/s)
kW      = 2.0;                       % turning gain
obsDist = 1.5;                       % LiDAR hits closer than this push the robot away (m)
kRep    = 0.3;                       % strength of that push

% ---- Occupancy grid (log-odds) ----
lFree = -0.4;                        % evidence added to cells a beam passes through
lOcc  =  0.85;                       % evidence added to the cell where a beam ends
lMax  =  5;                          % clamp so the map can still change its mind

% ---- Gas ----
gasSrc     = [17.5; 17];             % HIDDEN leak position (keep it outside obstacles)
gasPeak    = 100;                    % concentration at the source (ppm)
gasSigma   = 9;                      % how far the gas spreads (m)
gasNoise   = 0.3;                    % sensor noise (ppm, standard deviation)
sensorOff  = 0.4;                    % distance of front/back/left/right sensors from centre (m)
gradSmooth = 0.9;                    % 0..1, higher = smoother (slower) gradient estimate
gasFiltK   = 0.8;                    % 0..1, smoothing of the centre-sensor reading
foundLevel = 0.99;                   % stop when smoothed reading >= this fraction of peak

% ---- Display ----
plotEvery = 3;                       % update the figure every N steps
bgCol    = [0.07 0.08 0.11];         % background
fgCol    = [0.90 0.92 0.95];         % text / axes
gridCol  = [0.18 0.19 0.25];         % cell lines
colFree  = [0.35 0.75 1.00];         % known free space (sky blue)
colUnk   = [0.11 0.11 0.14];         % unknown space (charcoal)
colOcc   = [1.00 0.20 0.50];         % obstacle (hot pink)
colWall  = [0.85 0.85 0.90];         % walls drawn on the gas map
colPath  = [1.00 0.90 0.20];         % trajectory (yellow)
colRobot = [0.30 1.00 0.60];         % robot (mint)
robotLen = 0.8;                      % robot body length (m)  -> along the heading
robotWid = 0.5;                      % robot body width  (m)

%%  Robot surroundings
% Each row of obs is one rectangle: [x  y  width  height] (meters), (x,y) = bottom-left corner.
obs = [ 3    14.5  4    1   ;
        5     7.5  1    5   ;
        7.5  13.5  4.5  1   ;
        8.5   9.5  3    1   ;
        4     3.5  4    1   ;
       10     3.5  5    1   ;
       13.5   7    1    7.5 ];

trueGrid = false(N, N);              % true = obstacle.  row = y index, column = x index
for i = 1:size(obs,1)
    c1 = floor(obs(i,1)/res) + 1;   c2 = ceil((obs(i,1)+obs(i,3))/res);
    r1 = floor(obs(i,2)/res) + 1;   r2 = ceil((obs(i,2)+obs(i,4))/res);
    trueGrid(r1:r2, c1:c2) = true;
end
trueGrid([1 end], :) = true;         % outer walls
trueGrid(:, [1 end]) = true;

% The TRUE gas field. Gas spreads around walls, so concentration depends on the
% shortest path (geodesic distance) from the source, not the straight-line distance.
Gt = makeGasField(trueGrid, res, gasSrc, gasPeak, gasSigma);
xc = ((1:N) - 0.5)*res;              % x coordinate of each cell centre
yc = xc;

%% Initial State 
pose = start;                        % [x; y; heading]
traj = pose(1:2);                    % path history
L    = zeros(N, N);                  % log-odds occupancy map (0 = unknown)
P    = 0.5*ones(N, N);               % occupancy probability (0.5 = unknown)
Gsum = zeros(N, N);                  % gas map: sum of readings per cell
Gcnt = zeros(N, N);                  % gas map: number of readings per cell
gradEst = [0; 0];                    % smoothed gas gradient (world frame)
gasFilt = 0;                         % smoothed centre-sensor reading

readLog  = nan(1, maxSteps);         % logs for the metrics
distLog  = nan(1, maxSteps);
posLog   = nan(2, maxSteps);
minClear = inf;                      % smallest LiDAR distance seen during the run
found    = false;
nLog     = 0;

%%  Figure
fig = figure('Color', bgCol, 'Position', [30 80 1550 560]);
cmapOcc = interp1([0 0.5 1], [colFree; colUnk; colOcc], linspace(0,1,256));
cmapGas = turbo(256);

% --- Panel 1: true environment (what we know as the designer, not the robot) ---
ax1 = subplot(1,3,1);
imagesc(ax1, [0 worldSize], [0 worldSize], double(trueGrid));
colormap(ax1, [colUnk; colOcc]);
axis(ax1, 'xy', 'equal', 'tight');  hold(ax1, 'on');
drawGrid(ax1, worldSize, res, gridCol);
contour(ax1, xc, yc, Gt, [5 10 20 40 60 80 95], 'LineColor', [0.45 0.55 0.65], 'LineWidth', 0.8);
plot(ax1, gasSrc(1), gasSrc(2), 'x', 'Color', [1 0.3 0.3], 'MarkerSize', 12, 'LineWidth', 2);
styleAxes(ax1, bgCol, fgCol);
title(ax1, 'TRUE WORLD  (x = hidden leak, lines = true gas levels)', 'Color', fgCol);
hTraj1 = plot(ax1, traj(1,:), traj(2,:), '-', 'Color', colPath, 'LineWidth', 1.8);
hHits  = plot(ax1, NaN, NaN, '.', 'Color', [1 1 1], 'MarkerSize', 6);
[rx, ry] = robotShape(pose, robotLen, robotWid);
hRobot = patch(ax1, rx, ry, colRobot, 'EdgeColor', 'w', 'LineWidth', 1.2);   % rectangle robot
hHead  = plot(ax1, NaN, NaN, '-', 'Color', colRobot, 'LineWidth', 2);

% --- Panel 2: occupancy map built by the robot ---
ax2 = subplot(1,3,2);
hMap = imagesc(ax2, [0 worldSize], [0 worldSize], P);
colormap(ax2, cmapOcc);  caxis(ax2, [0 1]);
axis(ax2, 'xy', 'equal', 'tight');  hold(ax2, 'on');
drawGrid(ax2, worldSize, res, gridCol);
styleAxes(ax2, bgCol, fgCol);
title(ax2, 'LIDAR OCCUPANCY MAP', 'Color', fgCol);
hTraj2 = plot(ax2, traj(1,:), traj(2,:), '-', 'Color', colPath, 'LineWidth', 1.8);

% --- Panel 3: gas map built by the robot ---
ax3 = subplot(1,3,3);
hGas = imagesc(ax3, [0 worldSize], [0 worldSize], ...
               renderGasMap(Gsum, Gcnt, L, cmapGas, gasPeak, bgCol, colWall));
colormap(ax3, cmapGas);  caxis(ax3, [0 gasPeak]);
axis(ax3, 'xy', 'equal', 'tight');  hold(ax3, 'on');
drawGrid(ax3, worldSize, res, gridCol);
styleAxes(ax3, bgCol, fgCol);
title(ax3, 'GAS MAP', 'Color', fgCol);
hTraj3  = plot(ax3, traj(1,:), traj(2,:), '-', 'Color', [1 1 1], 'LineWidth', 1);
hRobot3 = patch(ax3, rx, ry, colRobot, 'EdgeColor', 'w', 'LineWidth', 1.2);
hEst    = plot(ax3, NaN, NaN, 'c+', 'MarkerSize', 14, 'LineWidth', 2);
cb = colorbar(ax3);  cb.Color = fgCol;
cb.Label.String = 'Gas concentration (ppm)';  cb.Label.Color = fgCol;

sgtitle('Autonomous Gas-Leak Localization Robot with LiDAR Occupancy Mapping', ...
        'Color', fgCol, 'FontWeight', 'bold');

%%  Main Loop: sense -> map -> think -> move -> draw 
for step = 1:maxSteps

    % ---------- SENSE ----------
    ranges = lidarScan(pose, trueGrid, res, relAng, maxR);          % 72 distances
    minClear = min(minClear, min(ranges));
    [gas, gasPts] = readGasSensors(pose, Gt, res, sensorOff, gasNoise);   % 5 readings

    % ---------- MAP ----------
    L = updateMap(L, pose, ranges, relAng, maxR, res, lFree, lOcc, lMax);
    P = 1 - 1./(1 + exp(L));                                         % log-odds -> probability
    [Gsum, Gcnt] = updateGasMap(Gsum, Gcnt, gasPts, gas, res);

    % ---------- ESTIMATE THE GAS GRADIENT ----------
    % gas = [centre front back left right]
    gradRobot = [gas(2) - gas(3);                 % front - back  -> gradient along robot's forward axis
                 gas(4) - gas(5)] / (2*sensorOff);% left  - right -> gradient along robot's left axis
    th = pose(3);
    Rot = [cos(th) -sin(th); sin(th) cos(th)];    % rotate from robot frame to world frame
    gradWorld = Rot * gradRobot;
    if step == 1
        gradEst = gradWorld;
        gasFilt = gas(1);
    else
        gradEst = gradSmooth*gradEst + (1 - gradSmooth)*gradWorld;  % low-pass filter (removes noise)
        gasFilt = gasFiltK*gasFilt   + (1 - gasFiltK)*gas(1);
    end

    % ---------- LOG ----------
    nLog = step;
    readLog(step) = gasFilt;
    distLog(step) = norm(pose(1:2) - gasSrc);
    posLog(:, step) = pose(1:2);

    % ---------- THINK + MOVE ----------
    if gasFilt >= foundLevel*gasPeak
        found = true;
    end
    if ~found
        [v, w] = gasController(pose, gradEst, ranges, relAng, obsDist, kRep, kW, vMax);
        pose(1) = pose(1) + v*cos(pose(3))*dt;                       % differential-drive kinematics
        pose(2) = pose(2) + v*sin(pose(3))*dt;
        pose(3) = atan2(sin(pose(3) + w*dt), cos(pose(3) + w*dt));
        traj(:, end+1) = pose(1:2); %#ok<SAGROW>
    end

    % ---------- DRAW ----------
    if mod(step, plotEvery) == 0 || found
        hitMask = ranges < maxR;
        a = pose(3) + relAng(hitMask);
        set(hHits, 'XData', pose(1) + ranges(hitMask).*cos(a), ...
                   'YData', pose(2) + ranges(hitMask).*sin(a));
        [rx, ry] = robotShape(pose, robotLen, robotWid);      % corners of the rotated rectangle
        set(hRobot,  'XData', rx, 'YData', ry);
        set(hRobot3, 'XData', rx, 'YData', ry);
        set(hHead, 'XData', [pose(1), pose(1) + 0.8*cos(pose(3))], ...
                   'YData', [pose(2), pose(2) + 0.8*sin(pose(3))]);
        set(hTraj1, 'XData', traj(1,:), 'YData', traj(2,:));
        set(hTraj2, 'XData', traj(1,:), 'YData', traj(2,:));
        set(hTraj3, 'XData', traj(1,:), 'YData', traj(2,:));
        set(hMap, 'CData', P);
        set(hGas, 'CData', renderGasMap(Gsum, Gcnt, L, cmapGas, gasPeak, bgCol, colWall));
        [~, iBest] = max(readLog(1:nLog));          % where the reading was highest so far
        set(hEst, 'XData', posLog(1, iBest), 'YData', posLog(2, iBest));
        ax2.Title.String = sprintf('LIDAR OCCUPANCY MAP   |   step %d', step);
        ax3.Title.String = sprintf('GAS MAP   |   reading %.1f ppm', gasFilt);
        drawnow limitrate;
    end

    if found
        break;
    end

    % Stuck detection: moved less than 0.3 m during the last 150 steps
    if step > 300 && norm(traj(:, end) - traj(:, end-150)) < 0.3
        warning('Robot seems stuck. Try changing kRep, obsDist or the leak position.');
        break;
    end
end

%%  Metrices
[~, iBest] = max(readLog(1:nLog));
estSrc  = posLog(:, iBest);                              % estimated leak position
locErr  = norm(estSrc - gasSrc);                         % localisation error (m)
pathLen = sum(sqrt(sum(diff(traj, 1, 2).^2, 1)));        % travelled distance (m)

known    = (L ~= 0);                                     % cells the robot has seen
estOcc   = (L > 0);                                      % robot believes "occupied"
explored = 100*nnz(known)/numel(L);                      % % of world seen
mapAcc   = 100*nnz(known & (estOcc == trueGrid))/max(nnz(known), 1);
tp  = nnz(known &  estOcc &  trueGrid);
fp  = nnz(known &  estOcc & ~trueGrid);
fn  = nnz(known & ~estOcc &  trueGrid);
mapIoU = tp/max(tp + fp + fn, 1);                        % obstacle IoU (0..1)

results = struct('found', found, 'steps', nLog, 'time_s', nLog*dt, 'path_m', pathLen, ...
                 'locError_m', locErr, 'minClearance_m', minClear, ...
                 'explored_pct', explored, 'mapAccuracy_pct', mapAcc, 'mapIoU', mapIoU);

% reveal the truth on the gas map
plot(ax3, gasSrc(1), gasSrc(2), 'x', 'Color', [1 1 1], 'MarkerSize', 14, 'LineWidth', 2);
set(hEst, 'XData', estSrc(1), 'YData', estSrc(2));
text(ax3, 0.5, 0.8, 'white x = true leak,  cyan + = robot estimate', ...
     'Color', fgCol, 'FontSize', 9);

fprintf('\n=============== GasLeakBot RESULTS ===============\n');
if found, fprintf('Leak found            : YES\n'); else, fprintf('Leak found            : NO\n'); end
fprintf('Time / steps          : %.1f s / %d steps\n', nLog*dt, nLog);
fprintf('Path length           : %.1f m\n', pathLen);
fprintf('True leak             : (%.1f, %.1f)\n', gasSrc(1), gasSrc(2));
fprintf('Estimated leak        : (%.1f, %.1f)\n', estSrc(1), estSrc(2));
fprintf('Localisation error    : %.2f m\n', locErr);
fprintf('Minimum clearance     : %.2f m (smallest LiDAR distance)\n', minClear);
fprintf('Area explored         : %.1f %%\n', explored);
fprintf('Map accuracy          : %.1f %% (on explored cells)\n', mapAcc);
fprintf('Obstacle IoU          : %.2f\n', mapIoU);
fprintf('==================================================\n');

% Second figure: gas reading and distance to the leak over time
t = (1:nLog)*dt;
fig2 = figure('Color', bgCol, 'Position', [100 100 950 380]);
axA = subplot(1,2,1);
plot(axA, t, readLog(1:nLog), 'Color', colPath, 'LineWidth', 1.6);  hold(axA, 'on');
yline(axA, foundLevel*gasPeak, '--', 'Color', colOcc);
styleAxes(axA, bgCol, fgCol);
xlabel(axA, 'Time (s)', 'Color', fgCol);  ylabel(axA, 'Gas reading (ppm)', 'Color', fgCol);
title(axA, 'Gas reading vs time', 'Color', fgCol);
axB = subplot(1,2,2);
plot(axB, t, distLog(1:nLog), 'Color', colFree, 'LineWidth', 1.6);
styleAxes(axB, bgCol, fgCol);
xlabel(axB, 'Time (s)', 'Color', fgCol);  ylabel(axB, 'Distance to true leak (m)', 'Color', fgCol);
title(axB, 'Distance to leak vs time', 'Color', fgCol);

%% Local Functions

function Gt = makeGasField(trueGrid, res, src, peak, sigma)
% True gas concentration in every cell.
% Step 1: geodesic distance D from the source, spreading only through free cells.
% Step 2: Gaussian falloff of that distance.
    [nR, nC] = size(trueGrid);
    D = inf(nR, nC);
    c = floor(src(1)/res) + 1;
    r = floor(src(2)/res) + 1;
    D(r, c) = 0;
    changed = true;
    while changed                                   % repeat until nothing improves
        Dold = D;
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0, continue; end
                stepCost = res*hypot(dr, dc);      % straight = res, diagonal = 1.41*res
                D = min(D, shiftInf(D, dr, dc) + stepCost);
            end
        end
        D(trueGrid) = inf;                          % gas cannot enter walls
        changed = any(D(:) ~= Dold(:));
    end
    Gt = peak*exp(-D.^2/(2*sigma^2));               % exp(-inf) = 0 inside walls
end

function B = shiftInf(A, dr, dc)
% B(r,c) = A(r+dr, c+dc); cells that fall outside the matrix become Inf.
    [nR, nC] = size(A);
    B  = inf(nR, nC);
    rs = max(1, 1-dr):min(nR, nR-dr);
    cs = max(1, 1-dc):min(nC, nC-dc);
    B(rs, cs) = A(rs + dr, cs + dc);
end

function ranges = lidarScan(pose, trueGrid, res, relAng, maxR)
% Ray marching: walk along each beam in 5 cm steps until it enters an occupied cell.
    [nR, nC] = size(trueGrid);
    r = 0:0.05:maxR;
    ranges = maxR*ones(1, numel(relAng));
    for k = 1:numel(relAng)
        a  = pose(3) + relAng(k);
        c  = floor((pose(1) + r*cos(a))/res) + 1;
        rw = floor((pose(2) + r*sin(a))/res) + 1;
        outside = c < 1 | c > nC | rw < 1 | rw > nR;
        c  = min(max(c, 1), nC);
        rw = min(max(rw, 1), nR);
        hit = outside | trueGrid(sub2ind([nR nC], rw, c));
        i = find(hit, 1);
        if ~isempty(i)
            ranges(k) = r(i);
        end
    end
end

function L = updateMap(L, pose, ranges, relAng, maxR, res, lFree, lOcc, lMax)
% Log-odds occupancy update: cells along a beam are free, the end cell is occupied.
    [nR, nC] = size(L);
    for k = 1:numel(relAng)
        a = pose(3) + relAng(k);
        r = 0:res:(ranges(k) - res);
        if ~isempty(r)
            c  = min(max(floor((pose(1) + r*cos(a))/res) + 1, 1), nC);
            rw = min(max(floor((pose(2) + r*sin(a))/res) + 1, 1), nR);
            idx = unique(sub2ind([nR nC], rw, c));
            L(idx) = L(idx) + lFree;
        end
        if ranges(k) < maxR
            c  = min(max(floor((pose(1) + ranges(k)*cos(a))/res) + 1, 1), nC);
            rw = min(max(floor((pose(2) + ranges(k)*sin(a))/res) + 1, 1), nR);
            L(rw, c) = L(rw, c) + lOcc;
        end
    end
    L = max(min(L, lMax), -lMax);
end

function [gas, pts] = readGasSensors(pose, Gt, res, off, noise)
% Five gas sensors: centre, front, back, left, right. Returns readings (ppm) and positions.
    th   = pose(3);
    fwd  = [cos(th); sin(th)];
    left = [-sin(th); cos(th)];
    c    = pose(1:2);
    pts  = [c, c + off*fwd, c - off*fwd, c + off*left, c - off*left];   % 2 x 5
    cx   = pts(1,:)/res + 0.5;                      % position in "grid index" units
    cy   = pts(2,:)/res + 0.5;
    gas  = interp2(Gt, cx, cy, 'linear', 0) + noise*randn(1, 5);
    gas  = max(gas, 0);                             % a sensor cannot read negative ppm
end

function [Gsum, Gcnt] = updateGasMap(Gsum, Gcnt, pts, gas, res)
% Store every reading in the grid cell where it was measured (average later).
    [nR, nC] = size(Gsum);
    for i = 1:numel(gas)
        c = min(max(floor(pts(1,i)/res) + 1, 1), nC);
        r = min(max(floor(pts(2,i)/res) + 1, 1), nR);
        Gsum(r, c) = Gsum(r, c) + gas(i);
        Gcnt(r, c) = Gcnt(r, c) + 1;
    end
end

function [v, w] = gasController(pose, gradEst, ranges, relAng, obsDist, kRep, kW, vMax)
% Potential field: the GAS GRADIENT attracts, nearby LiDAR hits repel.
    att = gradEst/(norm(gradEst) + 1e-9);           % unit vector toward stronger gas

    rep = [0; 0];
    for k = 1:numel(ranges)
        if ranges(k) < obsDist
            a   = pose(3) + relAng(k);
            rep = rep - (obsDist - ranges(k))/obsDist*[cos(a); sin(a)];
        end
    end

    total   = att + kRep*rep;
    desired = atan2(total(2), total(1));            % direction we want to go
    err     = atan2(sin(desired - pose(3)), cos(desired - pose(3)));   % heading error in (-pi, pi]

    front = min(ranges(abs(relAng) < pi/3));        % slow down if something is close ahead
    slow  = min(1, max(0.15, (front - 0.4)/1.0));

    w = kW*err;                                     % turn toward the desired direction
    v = vMax*max(0, cos(err))*slow;                 % go forward only if roughly facing it
end

function rgb = renderGasMap(Gsum, Gcnt, L, cmap, peak, bgCol, wallCol)
% Convert the gas map into an RGB picture: colour = average ppm, dark = not visited,
% light grey = walls the robot has detected with LiDAR.
    Gm  = Gsum./max(Gcnt, 1);
    nC  = size(cmap, 1);
    idx = min(max(round(Gm/peak*(nC - 1)) + 1, 1), nC);
    rgb = reshape(cmap(idx, :), [size(Gm) 3]);
    for ch = 1:3
        layer = rgb(:, :, ch);
        layer(Gcnt == 0) = bgCol(ch);
        layer(L > 0)     = wallCol(ch);
        rgb(:, :, ch) = layer;
    end
end

function [X, Y] = robotShape(pose, len, wid)
% Four corners of the rectangular robot body, rotated by the heading and moved to (x, y).
    body = [-len/2   len/2   len/2  -len/2;      % x of the 4 corners (robot frame)
            -wid/2  -wid/2   wid/2   wid/2];     % y of the 4 corners (robot frame)
    th  = pose(3);
    Rot = [cos(th) -sin(th); sin(th) cos(th)];   % rotation matrix
    pts = Rot*body + pose(1:2);                  % rotate, then translate
    X = pts(1, :);
    Y = pts(2, :);
end

function drawGrid(ax, worldSize, step, col)
% Thin cell lines (one plot call per direction, NaN separates the segments).
    t = 0:step:worldSize;
    X = [t; t; NaN(size(t))];
    Y = [zeros(size(t)); worldSize*ones(size(t)); NaN(size(t))];
    plot(ax, X(:), Y(:), '-', 'Color', col, 'LineWidth', 0.3);
    plot(ax, Y(:), X(:), '-', 'Color', col, 'LineWidth', 0.3);
end

function styleAxes(ax, bgCol, fgCol)
    set(ax, 'Color', bgCol, 'XColor', fgCol, 'YColor', fgCol, 'Box', 'on');
    xlabel(ax, 'X (m)', 'Color', fgCol);
    ylabel(ax, 'Y (m)', 'Color', fgCol);
end

