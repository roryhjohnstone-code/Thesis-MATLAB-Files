function frames = drawTurbine3D_advanced(ax, sensorStatus, animateFlag, exportFlag, frames, steps, strainData, ~, sensorSpanPos, faultMask)
% drawTurbine3D_advanced
% - Parametric blade mesh (span x chord)
% - Deflection by 1st flapwise mode ~ sin(pi*x/L) scaled by strain
% - Vertex colors from strain heatmap
% Inputs:
%   ax               : target axes
%   sensorStatus     : logical(5x1) which sensors flagged (for markers)
%   animateFlag      : bool
%   exportFlag       : bool
%   frames           : struct array for GIF recording (optional)
%   steps            : animation steps
%   strainData       : [numSensors x T]
%   time             : [1 x T]
%   sensorSpanPos    : [numSensors x 1], span locations in [0..1]
%   faultMask        : struct with fields .delamSpan=[a b], .crackSpan, .impactSpan (optional)

if nargin<9 || isempty(sensorSpanPos), sensorSpanPos = linspace(0.15,0.85,size(strainData,1)).'; end
if nargin<10, faultMask = struct; end
if isempty(ax) || ~isgraphics(ax), frames = []; return; end

% ---- Mesh params ----
L   = 3.0;          % blade half-length (m) (visual)
Nxs = 60;           % spanwise divisions
Nyc = 12;           % chordwise divisions
chordRoot = 0.45; chordTip = 0.15;
twistRoot = 0; twistTip = 18*pi/180;  % radians
xSpan = linspace(0, L, Nxs).';        % 0..L
chord = linspace(chordRoot, chordTip, Nxs).';
twist = linspace(twistRoot, twistTip, Nxs).';

% chord coordinate (−0.5..+0.5 scaled by local chord)
yc = linspace(-0.5,0.5,Nyc);
ycMat = repmat(yc, Nxs,1);
chordMat = repmat(chord,1,Nyc);
twistMat = repmat(twist,1,Nyc);
X = repmat(xSpan,1,Nyc);           % span axis
Y = chordMat .* ycMat .* cos(twistMat);
Z = chordMat .* ycMat .* sin(twistMat);

% ---- Compute a per-span strain line from sensors (interpolate)
% Use mean absolute strain at final time index (or choose index)
if isempty(strainData)
    spanStrain = zeros(Nxs,1);
else
    tIdx = size(strainData,2);  % latest
    sVal = mean(abs(strainData(:,max(1,tIdx-10):tIdx)),2,'omitnan'); % smooth a bit
    xi   = sensorSpanPos(:)*L;
    spanStrain = interp1(xi, sVal, xSpan, 'pchip', 'extrap');
    spanStrain = max(spanStrain, 0);
end
spanStrain = spanStrain / (max(spanStrain)+eps);  % 0..1 normalize

% ---- Mode-shape deflection (flapwise) scaled by strain
mode1 = sin(pi * X / L);        % Nx × Ny
deflAmp = 0.25;                  % scale (m per unit)
W = deflAmp * mode1 .* repmat(spanStrain,1,Nyc);  % flapwise out-of-plane
Z = Z + W;                       % deform surface

% ---- Color map by span strain
C = repmat(spanStrain,1,Nyc);   % vertex colors

% ---- Animation loop (rotation)
for t = 1:max(1,steps)
    cla(ax); hold(ax,'on'); axis(ax,'equal');
    view(ax, 35, 22); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'X (span)'); ylabel(ax,'Y (chord)'); zlabel(ax,'Z (flap)');
    xlim(ax, [0 L]); ylim(ax, [-0.6 0.6]*chordRoot); zlim(ax, [-0.4 0.6]);

    % Hub + tower (simple)
    plot3(ax,[0 0],[0 0],[-0.5 0.0],'k-','LineWidth',6);  % tower stump
    
    % Rotate about hub slowly
    theta = 2*pi*(t-1)/max(1,steps);
    Rz = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

    % Apply rotation to blade
    P = [X(:)'; Y(:)'; Z(:)'];
    Pr = Rz * P;  % 3×(Nxs*Nyc)
    Xr = reshape(Pr(1,:), Nxs, Nyc);
    Yr = reshape(Pr(2,:), Nxs, Nyc);
    Zr = reshape(Pr(3,:), Nxs, Nyc);

    % Draw blade with strain colors
    surf(ax, Xr, Yr, Zr, C, 'EdgeColor','none', 'FaceAlpha',0.95);
    colormap(ax, jet);
    clim(ax,[0 1]);
    cb = colorbar(ax,'Location', 'eastoutside');
    cb.Label.String ='Normalized Strain';


    % Sensor markers
    for k = 1:numel(sensorSpanPos)
        xs = sensorSpanPos(k)*L; ys = 0; zs = deflAmp*sin(pi*xs/L)*spanStrain(max(1,round(xs/L*Nxs)));
        p = Rz * [xs; ys; zs];
        mk = 'o'; clr = [0 0 0];
        if sensorStatus(k), mk = 's'; clr = [1 0 0]; end
        plot3(ax, p(1), p(2), p(3), mk, 'MarkerFaceColor', clr, 'MarkerEdgeColor','w','MarkerSize',8);
        text(ax, p(1), p(2), p(3)+0.03, sprintf('S%d',k),'Color','w','FontWeight','bold');
    end

    % Fault overlays (optional shaded span bands)
    if isfield(faultMask,'delamSpan')
        xs = faultMask.delamSpan*L;
        fill3(ax, [xs(1) xs(2) xs(2) xs(1)], [0 0 0 0], [0 0.02 0.02 0], [0 1 0], 'FaceAlpha',0.15, 'EdgeColor','none');
        text(ax, mean(xs), 0, 0.05, 'Delam','Color',[0 1 0]);
    end
    if isfield(faultMask,'crackSpan')
        xs = mean(faultMask.crackSpan)*L;
        plot3(ax, [xs xs], [-0.1 0.1], [0.05 0.05],'w-','LineWidth',3);
        text(ax, xs, 0.12, 0.05, 'Crack','Color','w');
    end
    if isfield(faultMask,'impactSpan')
        xs = mean(faultMask.impactSpan)*L;
        plot3(ax, xs, 0, 0.1, 'p', 'MarkerSize',14, 'MarkerFaceColor','y','MarkerEdgeColor','k');
        text(ax, xs, 0, 0.14, 'Impact','Color','y');
    end

    title(ax,'Wind Turbine Blade — Strain Heatmap');
    hold(ax,'off'); drawnow limitrate;

    if exportFlag
        frames(t) = getframe(ax); 
    end
    if ~animateFlag, break; end
end

