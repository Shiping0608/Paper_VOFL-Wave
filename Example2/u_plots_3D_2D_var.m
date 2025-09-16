function u_plots_3D_2D_var()


% -------------------------------------------------------------------------
%coscos;
for n = [2 4 6 8]
    data1 = load("./2D_Spatial_TS2_coscos/data_coscos_N2048_t"+num2str(n)+".mat");
    u3 = data1.u;
    xx = data1.xx;
    yy = data1.yy;
    figure
    surf(xx,yy,u3+1,u3)
    shading interp
    title("$t =\ $"+num2str(n),'Fontsize',18,'interpret','latex')
    hold on
    pcolor(xx,yy,u3)
    shading interp
    zticks([0 0.5 1 1.5])
    zticklabels({'-1','-0.5','0','0.5'})
    view(-45,23.5)
    xlim([-10 10])
    ylim([-10 10])
    if n == 2
        zlim([0 1.6]);
    else
        zlim([0 1.5]);
    end
    box on
    pause(0.1)
end