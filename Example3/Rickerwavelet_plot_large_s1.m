function Rickerwavelet_plot_large_s1()


% Q2v2 large s2
for t=0.1:0.1:0.5
    data = load(".\Q2v2_large_s1_dt5e_5" ...
     +"\data_RickerWavelet_N1000_t"+num2str(t)+".mat");
    N = length(data.xx(1,:));
    h = data.xx(1,2)-data.xx(1,1);
    x = [data.xx(1,:), data.xx(1,N)+h];
    y = [data.yy(:,1); data.yy(N,1)+h];
    U = zeros(N+1,N+1);
    U(1:N,1:N) = data.u;
    U(N+1,1:N+1) = U(1,:);
    U(1:N+1,N+1) = U(:,1);

    figure
    imagesc(x,y,abs(U));
    set(gca, 'YDir', 'normal')
    axis equal;
    title("$t =\,$"+num2str(t),'fontsize',24,'interpreter','latex')
    axis([0 2 0 2])
    %axis = caxis;
    caxis([1e-6 0.03])
    xticks([])
    yticks([])
    %yticks([0 0.5 1 1.5 2])
    %xlabel("$x$",'fontsize',20,'interpreter','latex')
    %ylabel("$t$",'fontsize',20,'interpreter','latex')
    %colormap gray
    numColors = 64;

    % Create a custom colormap from gray to white
    cmap = zeros(numColors, 3); % Initialize the colormap matrix
    cmap(:, 1) = linspace(0.5, 1, numColors); % Red component (from gray to white)
    cmap(:, 2) = linspace(0.5, 1, numColors); % Green component (from gray to white)
    cmap(:, 3) = linspace(0.5, 1, numColors); % Blue component (from gray to white)
    colormap(cmap);
%     if t == 0.5
%         hcb = colorbar;
%         hcb.Position(1) = hcb.Position(1)+0.05; 
%         hcb.Position(3) = 0.5*hcb.Position(3); % thiner colorbar
%     end

    max(max(U))
    pause(0.1)
end
