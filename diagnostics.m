function [grid] = diagnostics(rho,ux,uy,uz,grid)

%Compute E
gamma = sqrt(1+ux.^2+uy.^2+uz.^2);
KE = (gamma - 1).*rho;
fprintf("KE TOTAL: %1.12f\n",sum(sum(KE))*grid.dx/grid.E0);
grid.E_vs_t(grid.iter) = sum(sum(KE));

% Run only at select iterations:
if (mod ( grid.iter, grid.Output_interval ) == 0 || grid.iter == grid.NT)

    % Clear the figure
    clf()

    % Color
    %color = [ (grid.NT - grid.iter)/grid.NT,0,(grid.iter)/grid.NT];
    levels = 25;

    %Plot the diagnostic output comparison to fig1 GA Sod
    subplot(3,3,1)
    contourf(grid.x,grid.y,rho,levels)
    hold on
    title("Density")
    ylabel("y")
    xlabel("x")
    colorbar()

    subplot(3,3,2)
    contourf(grid.x,grid.y,ux,levels)
    hold on
    title("Momentum [x]")
    ylabel("y")
    xlabel("x")
    colorbar()

    subplot(3,3,3)
    contourf(grid.x,grid.y,KE,levels)
    title("Kinetic Energy")
    xlabel("x")
    ylabel("y")
    hold on
    colorbar()


    subplot(3,3,4)
    contourf(grid.x,grid.y,uy,levels)
    hold on
    title("Momentum [y]")
    ylabel("y")
    xlabel("x")
    colorbar()


    subplot(3,3,7)
    contourf(grid.x,grid.y,uz,levels)
    hold on
    title("Momentum [z]")
    ylabel("y")
    xlabel("x")
    colorbar()

    subplot(3,3,5)
    plot(grid.time_vec(1:grid.iter),grid.E_vs_t(1:grid.iter),"*")
    title("Kinetic Energy (t)")
    ylabel("KE")
    xlabel("t")

    subplot(3,3,6)
    vx = ux(1,:)./gamma(1,:);
    plot(grid.x,vx,"black")
    hold on
    vy = uy(1,:)./gamma(1,:);
    plot(grid.x,vy,"red")
     hold on
    vz = uz(1,:)./gamma(1,:);
    plot(grid.x,vz,"blue")
    title("Velocity (x = 0,y)")
    ylabel("v")
    xlabel("x")
    legend("vx","vy","vz")


    subplot(3,3,8)
    plot(grid.x,ux(1,:),"black")
    hold on
    plot(grid.x,uy(1,:),"red")
     hold on
    plot(grid.x,uz(1,:),"blue")
    title("Momentum (x = 0,y)")
    ylabel("v")
    xlabel("x")
    legend("ux","uy","uz")


    subplot(3,3,9)
    contourf(grid.x,grid.y,gamma,levels)
    title("gamma")
    ylabel("y")
    xlabel("x")
    colorbar()


    pause(0.01)

end
