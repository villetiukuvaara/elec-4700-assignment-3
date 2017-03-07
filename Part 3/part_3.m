%% Part 1: Electron Modelling
% This part is a modification of assignment 1. Here are some parameters for
% the simuation, including the width $W$ in the $y$-direction and the
% length $L$ in the $x$-direction.

clear all;
close all;

L = 200e-9;
W = 100e-9;
qe = -1.60217662e-19; % Charge on electron
Vx = 0.1; % Voltage applied across L (x = 0 is positive)
Vy = 0; % Voltage applied across W (y = 0 is positive)
density = 1e15*100^2; % Concentration of electrons in 1/m^2
m0 = 9.10938356e-31;
m = 0.26*m0;
T = 300;
k = 1.38064852e-23;
vth = sqrt(2*k*T/m);
l = vth*0.2e-12; % Mean free path
population_size = 30000;
plot_population = 10;
time_step = W/vth/100;
iterations = 200;
p_scat = 1 - exp(-time_step/0.2e-12);
v_pdf = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/m));
% Set to 1 to watch the movies,
% or to 0 to just see the final plots
show_movie = 0;

%%
% Here, the boundaries can be set to be specular or diffusive. If they
% are diffusive, the electrons bounce off at a random angle rather than
% one symmetrical about the normal with the boundary.
%
% The non-periodic top and bottom boundaries can be set to be either
% specular (1) or diffusive (0) with the following parameters:
top_specular = 0;
bottom_specular = 0;

%%
% When the given voltages are applied, the electric field components in the
% solid are (assuming that the fields are uniform):

Ex = Vx/L
Ey = Vy/W

%%
% The force on each electron is

Fx = qe*Ex
Fy = qe*Ey

%%
% For one time step, this increases the speed in each direction by

dvx = Fx*time_step/m;
dvy = Fy*time_step/m;
dvx = dvx.*ones(population_size,1);
dvy = dvy.*ones(population_size,1);

%%
% For the simulations, these arrays will hold information about the
% state of the system, including the positions, velocities, and
% temperatures.

% Each row corresponds to an electron with the positions and velocities
% [x y vx vy]
state = zeros(population_size, 4);
trajectories = zeros(iterations, plot_population*2);
temperature = zeros(iterations,1);
J = zeros(iterations,2); % Current density as [Jx Jy] rows

%%
% The relationship between electron drift current density and average
% carrier velocity is derived as follows. Let $v_x$ and $v_y$ be the
% velocity components in $x$ and $y$ for each of the $N$ particles in the simulation.
% The the average carrier velocties are $\bar{v_x} = 1/N \sum v_x$ and
% $\bar{v_y} = 1/N \sum v_y$. The electron concentration is $\rho =
% 10^{15}$ cm^2 and the charge is $q$. The electron drift current density
% components are
%
% $$ J_x = (q\rho) \left( \frac{1}{N} \right) \sum_{n=1}^N v_{x,n}
% = \frac{q\rho}{N} \sum_{n=1}^N v_{x,n}$$
%
% $$ J_y = (q\rho) \left( \frac{1}{N} \right) \sum_{n=1}^N v_{y,n}
% = \frac{q\rho}{N} \sum_{n=1}^N v_{y,n}$$
%
% These equations are used to plot the current density over time.

% Generate an initial population
for i = 1:population_size
    angle = rand*2*pi;
    state(i,:) = [L*rand W*rand random(v_pdf) random(v_pdf)];
end

figure(1);
subplot(3,1,1);
plot([],[]);
axis([0 L/1e-9 0 W/1e-9]);
title(sprintf('Trajectories for %d of %d Electrons (Part 3)',...
    plot_population, population_size));
xlabel('x (nm)');
ylabel('y (nm)');

figure(1);
subplot(3,1,2);
temperature_plot = animatedline;
title('Semiconductor Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');
grid on;

figure(1);
subplot(3,1,3);
current_plot = animatedline;
title('Drift Current Density J_x');
xlabel('Time (s)');
ylabel('Current density (A/m)');
grid on;

%%
% Run through the simulation:
for i = 1:iterations
    % Update the velocities
    
    state(:,3) = state(:,3) + dvx;
    state(:,4) = state(:,4) + dvy;
    
    %Update the positions
    state(:,1:2) = state(:,1:2) + time_step.*state(:,3:4);

    j = state(:,1) > L;
    state(j,1) = state(j,1) - L;
    
    j = state(:,1) < 0;
    state(j,1) = state(j,1) + L;
    
    j = state(:,2) > W;

    if(top_specular)
        state(j,2) = 2*W - state(j,2);
        state(j,4) = -state(j,4);
    else % Diffusive
        % The electron bounces off at a random angle
        state(j,2) = W;
        v = sqrt(state(j,3).^2 + state(j,4).^2);
        angle = rand([sum(j),1])*2*pi;
        state(j,3) = v.*cos(angle);
        state(j,4) = -abs(v.*sin(angle));
    end
    
    j = state(:,2) < 0;
    
    if(bottom_specular)
        state(j,2) = -state(j,2);
        state(j,4) = -state(j,4);
    else % Diffusive
        % The electron bounces off at a random angle
        state(j,2) = 0;
        v = sqrt(state(j,3).^2 + state(j,4).^2);
        angle = rand([sum(j),1])*2*pi;
        state(j,3) = v.*cos(angle);
        state(j,4) = abs(v.*sin(angle));
    end
    
    % Scatter particles
    j = rand(population_size, 1) < p_scat;
    state(j,3:4) = random(v_pdf, [sum(j),2]);
    
    % Record the temperature
    temperature(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/population_size;
    
    % Record positions and velocities for subset of particles that will be graphed
    for j=1:plot_population
        trajectories(i, (2*j):(2*j+1)) = state(j, 1:2);
    end
    
    % Calculate and record the current density using the formula from above
    J(i, 1) = qe.*density.*mean(state(:,3));
    J(i, 2) = qe.*density.*mean(state(:,4));

    % Plot the temperature and current
    addpoints(temperature_plot, time_step.*i, temperature(i));
    addpoints(current_plot, time_step.*i, J(i,1));

    if(show_movie && mod(i,5) == 0)
        figure(1);
        subplot(3,1,1);
        hold off;
        plot(state(1:plot_population,1)./1e-9, state(1:plot_population,2)./1e-9, 'o');
        axis([0 L/1e-9 0 W/1e-9]);
        hold on;
        title(sprintf('Trajectories for %d of %d Electrons (Part 3)',...
        plot_population, population_size));
        xlabel('x (nm)');
        ylabel('y (nm)');
        pause(0.05);
    end
end

% Show trajectories after the movie is over
figure(1);
subplot(3,1,1);
title(sprintf('Electron Trajectories for %d of %d Electrons (Part 3)',...
    plot_population, population_size));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 L/1e-9 0 W/1e-9]);
grid on;
hold on;
for i=1:plot_population
    plot(trajectories(:,i*2)./1e-9, trajectories(:,i*2+1)./1e-9, '.');
end

%%
% The plot of current density over time shows that it starts at 0 and increases
% over a period of about 0.5 ps to a steady state current. For this case
% with 0.1 V, the steady state current is about 87 A/mm. The reason that it
% starts at 0 is that the velocity is distributed using the equilibrium
% velocity distribution at the start, with none of the particles having
% accelerated yet.

%%
% First, plot the electron density. This involves performing some smoothing
% (filtering) using the conv2 function (Gaussian filtering). 

density = hist3(state(:,1:2),[200 100])';
N = 20;
sigma = 1.5;
[x y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(2);
density = conv2(density,f,'same');
density = density/(W./size(density,1)*L./size(density,2));
surf(conv2(density,f,'same'));
title('Electron Density');
xlabel('x (nm)');
ylabel('y (nm)');

%%
% The temperature map is created using a similar procudure. The electron
% velocities are put into bins over space to calculate the temperature at
% different points:
temp_sum_x = zeros(ceil(L/1e-9),ceil(W/1e-9));
temp_sum_y = zeros(ceil(L/1e-9),ceil(W/1e-9));
temp_num = zeros(ceil(L/1e-9),ceil(W/1e-9));

% Look at velocities of all the particles
for i=1:population_size
    % Find which "bin" it belongs in:
    x = floor(state(i,1)/1e-9);
    y = floor(state(i,2)/1e-9);
    if(x==0)
        x = 1;
    end
    if(y==0)
        y= 1;
    end
    
    % Add its velocity components to the cumulative count:
    temp_sum_y(x,y) = temp_sum_y(x,y) + state(i,3)^2;
    temp_sum_x(x,y) = temp_sum_x(x,y) + state(i,4)^2;
    temp_num(x,y) = temp_num(x,y) + 1;
end

%%
% Now, with the velocities added up, calculate the temperatures:
temp = (temp_sum_x + temp_sum_y).*m./k./2./temp_num;
temp(isnan(temp)) = 0;
temp = temp';

%%
% Like with the density map, perform some smoothing:

N = 20;
sigma = 1.5;
[x y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(3);
surf(conv2(temp,f,'same'));
title('Temperature Map');
xlabel('x (nm)');
ylabel('y (nm)');