function sol_abs_plot_L128var()
% SOL_ABS_PLOT_L128VAR   Plot solution magnitude.
%
%   This function loads precomputed solution data from a .mat file and
%   generates a log-scaled density plot of |u(x,t)| over space-time.
%   The data corresponds to the test case:
%       s(x) = 1 + 0.001 * [1 - sin(pi * x / 64)]^7
%
%   Input:
%       None.
%
%   Data:
%       ./TS2_0.001sin7_L128_h64/data_sin_N<N>_t40.mat
%       containing:
%           solution  – Solution matrix, size (Nt x Nx)
%           x         – Spatial grid points
%           T0, T1    – Start and end times
%
%   Output:
%       A mesh plot (in 2D view) of |u(x,t)| with log-scaled colors.
%
%   Author:
%       Shiping Zhou
%
%   Date last updated:
%       September 15, 2025 (Monday)
% -------------------------------------------------------------------------

N = 16384;

% -------------------------------------------------------------------------
data = load(".\TS2_0.001sin7_L128_h64\data_sin_N"+num2str(N)+"_t40.mat");

u = data.solution;
x = data.x;
T0 = data.T0;
T1 = data.T1;

dtt = 0.1;
t = T0:dtt:T1;

[xx,tt] = meshgrid(x,t);

figure
mesh(xx,tt,abs(u'));
xlim([-64,64]);
view(2)
ax = gca;
axpos = ax.Position;
colormap turbo
% customize colorbar
hcb = colorbar('eastoutside');
    
set(gca,'ColorScale','log')
title("$s(x)=1+0.001[1-\sin(\pi x/64)]^7$",'interpreter','latex')
xlabel("$x$",'fontsize',14,'interpreter','latex')
ylabel('$t$','fontsize',14,'interpreter','latex')
pause(0.1);
