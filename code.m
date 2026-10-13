
%% Stokes equations on a staggered MAC grid
% Manufactured solution and second-order convergence
clear; clc; close all;

Ns = [8 16 32 64 128];
errors = zeros(numel(Ns),3);

%% Solve and calculate errors
for k = 1:numel(Ns)
    N = Ns(k);

    [U,V,P,Uex,Vex,Pex] = solveStokes(N);

    errors(k,1) = norm(U(:)-Uex(:))/norm(Uex(:));
    errors(k,2) = norm(V(:)-Vex(:))/norm(Vex(:));
    errors(k,3) = norm(P(:)-Pex(:))/norm(Pex(:));
end

% Observed convergence orders
orders = log2(errors(1:end-1,:)./errors(2:end,:));

fprintf('Convergence orders [U, V, P]:\n');
disp(orders);


%% Figure 1: Velocity magnitude and streamlines

N = Ns(end);
[U,V] = solveStokes(N);

h = 1/N;
xc = ((1:N)-0.5)*h;

% Interpolate velocities to pressure centers
Uc = 0.5*(U+circshift(U,1,1));
Vfull = [-3.5*ones(N,1), V, -3.5*ones(N,1)];
Vc = 0.5*(Vfull(:,1:end-1)+Vfull(:,2:end));

Uc = Uc';
Vc = Vc';

% Periodic extension
xext = [xc(end)-1, xc, xc(1)+1];
yext = [0, xc, 1];

Uext = [Uc(:,end), Uc, Uc(:,1)];
Vext = [Vc(:,end), Vc, Vc(:,1)];

Uext = [zeros(1,N+2); Uext; zeros(1,N+2)];
Vext = [-3.5*ones(1,N+2); Vext; ...
        -3.5*ones(1,N+2)];

% Fine visualization grid
xq = linspace(0,1,301);
yq = linspace(0,1,301);
[Xq,Yq] = meshgrid(xq,yq);

Uq = interp2(xext,yext,Uext,Xq,Yq,'linear');
Vq = interp2(xext,yext,Vext,Xq,Yq,'linear');

speed = hypot(Uq,Vq);

%% Plot
fig1 = figure('Color','w','Position',[100 100 760 670]);

contourf(Xq,Yq,speed,100,'LineStyle','none');
hold on;

% Colorblind-safe cividis-style palette
anchors = [
    0.000 0.135 0.305
    0.170 0.286 0.430
    0.340 0.410 0.485
    0.525 0.540 0.465
    0.735 0.680 0.385
    0.995 0.910 0.220
];

cmap = interp1(linspace(0,1,size(anchors,1)), ...
    anchors,linspace(0,1,256));

colormap(fig1,cmap);
clim([1.5 5.5]);

cb = colorbar;
cb.Color = [0.1 0.1 0.1];
cb.FontSize = 12;
cb.Ticks = 1.5:0.5:5.5;

%% Streamlines
sx = linspace(0.025,0.975,25);
sy = 0.995*ones(size(sx));

streams = stream2(Xq,Yq,Uq,Vq,sx,sy);

% Draw streamlines
hl = streamline(streams);
set(hl,'Color','w','LineWidth',1.3);

% Add directional arrows along each streamline
for k = 1:numel(streams)
    pts = streams{k};

    if size(pts,1) < 10
        continue;
    end

    % Arrows at approximately 30% and 70% of each path
    idx = round([0.3 0.7]*size(pts,1));

    for j = idx
        x0 = pts(j,1);
        y0 = pts(j,2);

        % Interpolate local velocity direction
        ux = interp2(Xq,Yq,Uq,x0,y0);
        vy = interp2(Xq,Yq,Vq,x0,y0);

        mag = hypot(ux,vy);

        if mag > 0
            % Arrow length in coordinate units
            len = 0.035;

            quiver(x0,y0, ...
                len*ux/mag,len*vy/mag,0, ...
                'Color','w', ...
                'LineWidth',1.2, ...
                'MaxHeadSize',1.8, ...
                'AutoScale','off');
        end
    end
end

%% Formatting
axis equal;
axis([0 1 0 1]);

xticks(0:0.2:1);
yticks(0:0.2:1);

xlabel('$x$','Interpreter','latex','FontSize',17, ...
    'Color','k');
ylabel('$y$','Interpreter','latex','FontSize',17, ...
    'Color','k');
title('$\|\mathbf{u}\|$','Interpreter','latex', ...
    'FontSize',20,'Color','k');

set(gca,'Color','w', ...
    'XColor','k','YColor','k', ...
    'FontSize',13,'LineWidth',0.9, ...
    'Box','on','Layer','top');

% Save PDFs next to this script, not in MATLAB's current folder
outdir = fileparts(mfilename('fullpath'));
if isempty(outdir), outdir = pwd; end

exportgraphics(fig1, fullfile(outdir,'stokes_velocity.pdf'), ...
    'ContentType','vector', 'BackgroundColor','white');



%% Figure 2: Error convergence

% Okabe-Ito colorblind-safe colors
blue   = [0.000 0.447 0.698];
green  = [0.000 0.620 0.451];
orange = [0.902 0.624 0.000];
black  = [0.12 0.12 0.12];

fig2 = figure('Color','w','Position',[150 150 850 580]);

ax = axes(fig2);
ax.XScale = 'log';
ax.YScale = 'log';
hold(ax,'on');

% Explicitly set light background and dark axes
set(ax,'Color','w', ...
    'XColor','k','YColor','k', ...
    'FontSize',13, ...
    'LineWidth',0.9, ...
    'Box','on');

%% Numerical errors
loglog(Ns,errors(:,1),'-o', ...
    'Color',blue,'LineWidth',1.6, ...
    'MarkerSize',7,'MarkerFaceColor',blue);

loglog(Ns,errors(:,2),'-^', ...
    'Color',green,'LineWidth',1.6, ...
    'MarkerSize',7,'MarkerFaceColor',green);

loglog(Ns,errors(:,3),'-p', ...
    'Color',orange,'LineWidth',1.6, ...
    'MarkerSize',8,'MarkerFaceColor',orange);

%% Reference convergence slopes
Nref = linspace(min(Ns),max(Ns),300);

% First order: O(N^-1)
ref1 = 1.5*errors(1,1)*(Ns(1)./Nref);

% Second order: O(N^-2)
ref2 = 0.7*errors(1,1)*(Ns(1)./Nref).^2;

loglog(Nref,ref1,'k--','LineWidth',1.5);
loglog(Nref,ref2,'-', ...
    'Color',black,'LineWidth',1.5);

%% Labels and title
xlabel('$N$','Interpreter','latex', ...
    'FontSize',16,'Color','k');

ylabel('Relative $L^2$ Error', ...
    'Interpreter','latex', ...
    'FontSize',16,'Color','k');

title('Error Plots for $U$, $V$, and $P$', ...
    'Interpreter','latex', ...
    'FontSize',18,'Color','k');

%% Legend
lgd = legend({ ...
    'Error $U$', ...
    'Error $V$', ...
    'Error $P$', ...
    'First Order Check', ...
    'Second Order Check'}, ...
    'Interpreter','latex', ...
    'Location','southwest', ...
    'FontSize',12);

lgd.Color = 'w';
lgd.TextColor = 'k';
lgd.EdgeColor = [0.3 0.3 0.3];

%% Grid and axes
grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

ax.GridColor = [0.6 0.6 0.6];
ax.GridAlpha = 0.25;
ax.MinorGridAlpha = 0.12;

xlim([min(Ns) max(Ns)]);

% Export
exportgraphics(fig2, fullfile(outdir,'stokes_errors.pdf'), ...
    'ContentType','vector', 'BackgroundColor','white');






%% =====================================================
% MAC grid Stokes solver
% ======================================================
function [U,V,P,Uex,Vex,Pex] = solveStokes(N)

h = 1/N;
I = speye(N);
M = N-1;

% MAC grids (x is the first array dimension)
xc = ((1:N)-0.5)*h;

[Xu,Yu] = ndgrid((1:N)*h,xc);
[Xv,Yv] = ndgrid(xc,(1:M)*h);
[Xp,Yp] = ndgrid(xc,xc);

% Manufactured solutions
Uex = sin(2*pi*Yu).*sin(2*pi*Xu);

Vex = -3.5 + cos(2*pi*Xv).* ...
      (cos(2*pi*Yv)-1);

Pex = sin(2*pi*Yp).*cos(2*pi*Xp);

% Forcing functions
f = (-8*pi^2+2*pi)* ...
    sin(2*pi*Yu).*sin(2*pi*Xu);

g = cos(2*pi*Xv).* ...
    (4*pi^2-(8*pi^2+2*pi)*cos(2*pi*Yv));

%% Periodic horizontal operators
S = sparse(1:N,[2:N 1],1,N,N);

Txx = (S+S'-2*I)/h^2;
Gx1 = (S-I)/h;
Dx1 = -Gx1';

%% Vertical Laplacians
e = ones(N,1);

Tyu = spdiags([e -2*e e],-1:1,N,N)/h^2;
Tyu(1,1) = -3/h^2;
Tyu(end,end) = -3/h^2;

ev = ones(M,1);
Tyv = spdiags([ev -2*ev ev],-1:1,M,M)/h^2;

%% Two-dimensional operators
Lu = kron(I,Txx)+kron(Tyu,I);
Lv = kron(speye(M),Txx)+kron(Tyv,I);

Gx = kron(I,Gx1);
Dx = kron(I,Dx1);

Gy1 = sparse([1:M 1:M],[1:M 2:N], ...
             [-ones(1,M) ones(1,M)]/h,M,N);

Gy = kron(Gy1,I);
Dy = -Gy';

%% Boundary contributions
Vwall = -3.5;

gbc = zeros(N,M);
gbc(:,1) = gbc(:,1)-Vwall/h^2;
gbc(:,end) = gbc(:,end)-Vwall/h^2;

hbc = zeros(N,N);
hbc(:,1) = Vwall/h;
hbc(:,end) = -Vwall/h;

%% Assemble system
nu = N*N;
nv = N*M;
np = N*N;

A = [Lu, sparse(nu,nv), -Gx;
     sparse(nv,nu), Lv, -Gy;
     Dx, Dy, sparse(np,np)];

rhs = [f(:);g(:)+gbc(:);hbc(:)];

% Zero-mean pressure constraint
c = sparse([zeros(nu+nv,1);ones(np,1)]);

A = [A c; c' 0];
rhs = [rhs;0];

%% Solve
sol = A\rhs;

U = reshape(sol(1:nu),N,N);
V = reshape(sol(nu+1:nu+nv),N,M);
P = reshape(sol(nu+nv+1:nu+nv+np),N,N);

% Remove pressure constant
P = P-mean(P(:));
Pex = Pex-mean(Pex(:));

end
