
%% Porous microstructure: MAC-grid Brinkman simulation
clear; clc; close all;

N = 96;
h = 1/N;
mu = 1;
alphaSolid = 1e5;

% Applied force: downward
Fx = 0;
Fy = -1;

% Circular grains: [x-center, y-center, radius]
grains = [
    0.20 0.20 0.105
    0.52 0.17 0.115
    0.83 0.23 0.095
    0.34 0.49 0.125
    0.70 0.51 0.110
    0.14 0.78 0.100
    0.48 0.82 0.115
    0.83 0.81 0.110
];

%% MAC-grid coordinates (periodic in x and y)
xc = ((1:N)-0.5)*h;
xf = (0:N-1)*h;

[Xu,Yu] = ndgrid(xf,xc);  % U on vertical faces
[Xv,Yv] = ndgrid(xc,xf);  % V on horizontal faces
[Xp,Yp] = ndgrid(xc,xc);  % P at cell centers

%% Geometry masks on all three grids
solidU = grainMask(Xu,Yu,grains);
solidV = grainMask(Xv,Yv,grains);
solidP = grainMask(Xp,Yp,grains);

alphaU = alphaSolid*double(solidU);
alphaV = alphaSolid*double(solidV);

porosity = 1-mean(solidP(:));

%% Periodic one-dimensional operators
I = speye(N);
S = sparse(1:N,[2:N 1],1,N,N);

T = (S+S'-2*I)/h^2;

% Gradient from cell centers to faces
G1 = (I-S')/h;

% Divergence from faces to cell centers
D1 = -G1';

%% Two-dimensional Laplacian
L = kron(I,T)+kron(T,I);

%% Gradient and divergence
Gx = kron(I,G1);
Gy = kron(G1,I);

Dx = -Gx';
Dy = -Gy';

%% Brinkman velocity operators
nu = N^2;
np = N^2;

Lu = mu*L-spdiags(alphaU(:),0,nu,nu);
Lv = mu*L-spdiags(alphaV(:),0,nu,nu);

%% Assemble Stokes-Brinkman system
A = [
    Lu, sparse(nu,nu), -Gx;
    sparse(nu,nu), Lv, -Gy;
    Dx, Dy, sparse(np,np)
];

rhs = [
    -Fx*ones(nu,1);
    -Fy*ones(nu,1);
    zeros(np,1)
];

% Pressure has an arbitrary additive constant
c = sparse([zeros(2*nu,1);ones(np,1)]);
A = [A,c;c',0];
rhs = [rhs;0];

%% Solve
fprintf('Solving %d unknowns...\n',size(A,1));
tic;
sol = A\rhs;
elapsed = toc;
fprintf('Solve time: %.2f seconds\n',elapsed);

U = reshape(sol(1:nu),N,N);
V = reshape(sol(nu+1:2*nu),N,N);
P = reshape(sol(2*nu+1:3*nu),N,N);
P = P-mean(P(:));

%% Check numerical residuals
ru = Lu*U(:)-Gx*P(:)+Fx;
rv = Lv*V(:)-Gy*P(:)+Fy;
rd = Dx*U(:)+Dy*V(:);

fprintf('Momentum residual U: %.2e\n',norm(ru,inf));
fprintf('Momentum residual V: %.2e\n',norm(rv,inf));
fprintf('Divergence residual: %.2e\n',norm(rd,inf));

%% Cell-centered velocity for plotting
% Periodic averaging of neighboring face values
Uc = 0.5*(U+circshift(U,-1,1));
Vc = 0.5*(V+circshift(V,-1,2));

speed = hypot(Uc,Vc);

% Superficial velocity: average over total domain
Umean = mean(U(:));
Vmean = mean(V(:));

% Permeability in direction of applied force
Kyy = mu*Vmean/Fy;

fprintf('Porosity: %.4f\n',porosity);
fprintf('Mean velocity: (%.4e, %.4e)\n',Umean,Vmean);
fprintf('Estimated Kyy: %.6e\n',Kyy);

%% Figure 1: microstructure
figure('Color','w');
imagesc(xc,xc,solidP');
axis xy equal tight;
colormap(gca,flipud(gray));
xlabel('x'); ylabel('y');
title(sprintf('Porous microstructure, porosity = %.3f',porosity));


%% Figure 2: velocity magnitude with directional streamlines

% Cell-centered velocities (from the MAC grid)
Uc = 0.5*(U + circshift(U,-1,1));
Vc = 0.5*(V + circshift(V,-1,2));

% Interpolate onto a finer grid for smooth visualization
nFine = 450;
xq = linspace(0,1,nFine);
yq = linspace(0,1,nFine);
[Xq,Yq] = meshgrid(xq,yq);

% Periodic extension of cell-centered data
xext = [xc(end)-1, xc, xc(1)+1];

Uext = [Uc(end,:); Uc; Uc(1,:)];
Vext = [Vc(end,:); Vc; Vc(1,:)];

Uext = [Uext(:,end), Uext, Uext(:,1)];
Vext = [Vext(:,end), Vext, Vext(:,1)];

yext = [xc(end)-1, xc, xc(1)+1];

[Xext,Yext] = ndgrid(xext,yext);

Uq = interp2(Xext',Yext',Uext',Xq,Yq,'linear');
Vq = interp2(Xext',Yext',Vext',Xq,Yq,'linear');

speedq = hypot(Uq,Vq);

% Solid mask on the fine grid
solidFine = grainMask(Xq,Yq,grains);

% Prevent streamlines from entering obstacles
Ustream = Uq;
Vstream = Vq;
Ustream(solidFine) = NaN;
Vstream(solidFine) = NaN;

speedq(solidFine) = NaN;

%% Colorblind-safe cividis-inspired colormap
anchors = [
    0.000 0.135 0.304
    0.135 0.258 0.427
    0.275 0.381 0.465
    0.408 0.490 0.459
    0.565 0.600 0.425
    0.745 0.714 0.350
    0.900 0.820 0.240
    0.996 0.909 0.218
];

t0 = linspace(0,1,size(anchors,1));
t1 = linspace(0,1,256);

cmap = interp1(t0,anchors,t1,'pchip');
cmap = min(max(cmap,0),1);

%% Create figure
figPorous = figure('Color','w', ...
    'Position',[100 100 1000 800]);

ax = axes(figPorous);
hold(ax,'on');

% Smooth velocity magnitude contours
contourf(ax,Xq,Yq,speedq,45, ...
    'LineStyle','none');

colormap(ax,cmap);

cb = colorbar(ax);
cb.Label.String = 'Velocity magnitude';
cb.FontSize = 12;
cb.Color = 'k';

% Use a fixed, meaningful color range
clim(ax,[0 max(speedq(:),[],'omitnan')]);

%% Draw obstacles
theta = linspace(0,2*pi,300);

for k = 1:size(grains,1)
    cx = grains(k,1);
    cy = grains(k,2);
    r  = grains(k,3);

    fill(ax,cx+r*cos(theta),cy+r*sin(theta), ...
        [0.19 0.23 0.27], ...
        'EdgeColor',[0.08 0.10 0.13], ...
        'LineWidth',1.2);
end

%% Compute streamlines
% Seeds are placed near the top of the domain
sx = linspace(0.02,0.98,36);
sy = 0.97*ones(size(sx));

% Remove seeds inside obstacles
valid = ~grainMask(sx,sy,grains);
sx = sx(valid);
sy = sy(valid);

streams = stream2(Xq,Yq,Ustream,Vstream,sx,sy, ...
    [0.3 12000]);

%% Plot streamlines with directional arrows
for k = 1:numel(streams)

    path = streams{k};

    if isempty(path)
        continue;
    end

    % Keep only the finite portion of the streamline
    valid = all(isfinite(path),2);

    % Stop at the first invalid point to avoid
    % connecting disconnected streamline segments
    firstInvalid = find(~valid,1,'first');

    if ~isempty(firstInvalid)
        path = path(1:firstInvalid-1,:);
    end

    if size(path,1)<12
        continue;
    end

    % Remove consecutive duplicate points
    ds = hypot(diff(path(:,1)),diff(path(:,2)));
    keep = [true; ds > 1e-12];
    path = path(keep,:);

    if size(path,1)<2
        continue;
    end

    % Compute cumulative arc length
    ds = hypot(diff(path(:,1)),diff(path(:,2)));
    s = [0;cumsum(ds)];

    if ~isfinite(s(end)) || s(end)<0.06
        continue;
    end

    % Draw streamline
    plot(ax,path(:,1),path(:,2), ...
        'Color','w','LineWidth',1.1);

    % Arrow positions along streamline
    arrowLocations = 0.09:0.15:(s(end)-0.03);

    for j = 1:numel(arrowLocations)

        sj = arrowLocations(j);

        % Interpolate arrow position
        px = interp1(s,path(:,1),sj,'linear');
        py = interp1(s,path(:,2),sj,'linear');

        % Find local streamline direction
        idx = find(s<=sj,1,'last');

        if isempty(idx) || idx>=size(path,1)
            continue;
        end

        tx = path(idx+1,1)-path(idx,1);
        ty = path(idx+1,2)-path(idx,2);

        tn = hypot(tx,ty);

        if tn<1e-12
            continue;
        end

        tx = tx/tn;
        ty = ty/tn;

        % Draw directional arrow
        arrowLength = 0.025;

        quiver(ax,px,py, ...
            arrowLength*tx,arrowLength*ty,0, ...
            'Color','w', ...
            'LineWidth',1.1, ...
            'MaxHeadSize',2.5, ...
            'AutoScale','off');
    end
end

%% Formatting
axis(ax,'equal');
axis(ax,[0 1 0 1]);

ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.FontSize = 13;
ax.LineWidth = 1.1;
ax.Box = 'on';
ax.Layer = 'top';

xlabel(ax,'x','FontSize',15);
ylabel(ax,'y','FontSize',15);

title(ax,'Velocity magnitude and streamlines', ...
    'FontSize',16,'FontWeight','normal');

%% Export vector PDF
% Save PDF next to this script, not in a hardcoded path
outdir = fileparts(mfilename('fullpath'));
if isempty(outdir), outdir = pwd; end

exportgraphics(figPorous, fullfile(outdir,'porous_velocity.pdf'), ...
    'ContentType','vector', ...
    'BackgroundColor','white');


%% Local function: periodic circular-grain mask
function solid = grainMask(X,Y,grains)
    solid = false(size(X));

    for k = 1:size(grains,1)
        cx = grains(k,1);
        cy = grains(k,2);
        r = grains(k,3);

        % Minimum-image distances on periodic unit square
        dx = mod(X-cx+0.5,1)-0.5;
        dy = mod(Y-cy+0.5,1)-0.5;

        solid = solid | (dx.^2+dy.^2 <= r^2);
    end
end
