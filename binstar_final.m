clear all
close all

% Constants
G  = 6.67e-11;                      % Gravitational constant (N-m^2/kg^2)
earth_year = 365.25*24*60*60;       % Duration of year (s)
%% Circular Orbit -- Verified and correct
 binstar(6.0e24, 6.0e24, 10e8, 1.5e8,  50e6, 8*180000,    earth_year/4, 'b', 'r', 'm', 'mo')
 binstar(6.0e24, 6.0e24, 10e8, 1.5e8,  50e6, 16*180000,    earth_year/4, 'b', 'r', 'k', 'ko') %48e7
% binstar(6.0e24, 6.0e24, 10e8, 3e8, 9e7, 2*180000,  earth_year/4, 'b', 'r', 'k', 'ko' )
% hold off
%% Elliptical Orbit -- Planet apparently gets stripped by star2 
one_loop = earth_year/8 * .79
% binstar(6.0e24, 6.0e24, 10e8, 1.5e8, 9e7,     earth_year/60,  earth_year/8, 'b', 'r', 'k', 'ko' )
% hold off: 
% Planet gets stripped twice 
binstar(6.0e24, 6.0e24, 10e8, 1.5e8, 9e7,    4*earth_year/50,  8*one_loop, 'b', 'r', 'k', 'ko')
%% Planet doesn't get stripped twice 
binstar(6.0e24, 6.0e24, 10e8, 1.5e8, 90e5,    4*earth_year/50,  8*one_loop, 'b', 'r', 'k', 'ko')

function binstar(m1, m2, a,b,rp, clockmax, tmax, c1,c2,c3,c4)
    G  = 6.67e-11; 
    M     = m1 + m2;                    % Combined mass of stars
    va    = sqrt(b/a * (2*G*M)/(a+b));  % Initial velocity of stars
    U     = [0, va, 0];                 % Initial velocity vector    
    s1pos = -m2/M * [a, 0, 0];          % Position of Star1
    U1    = -m2/M * U;                  % Velocity of star 1
    s2pos =  m1/M * [a, 0, 0];          % Position of Star2
    U2    =  m1/M * U;                  % Velocity of star 2
    vpi   = sqrt(G*m1/rp);              % Initial velocity of planet around Star 1
    ppos  = s1pos + [rp, 0, 0];         % Initial position of planet
    UP    = U1 + [0, vpi, 0];           % Initial velocity of planet in system
    X = s2pos - s1pos;                  % X = X2-X1 | Distance btwn stars›
    com = (m1*s1pos + m2*s2pos)/(M);    % Center of mass
    
    % Create handles for stars, planet & trajectory:
    hc = plot(com(1), com(2), 'co', 'linewidth', 2);       % plot center of mass
    hold on
    title(['Circular Orbit:', newline, 'rp = ', num2str(rp,'%.2e'), ' a = ', num2str(a,'%.2e'), ' b =', num2str(b, '%.2e')],...
    'Color','black');
    hp1 = plot(s1pos(1), s1pos(2), 'b*','linewidth', 5);   % star1
    ht1 = plot(s1pos(1), s1pos(2),  c1 ,'linewidth', 2);   % star1 trajectory
    hp2 = plot(s2pos(1), s2pos(2), 'r*','linewidth', 5);   % star2
    ht2 = plot(s2pos(1), s2pos(2),  c2 ,'linewidth', 2);   % star2 trajectory
    hpp = plot(ppos(1), ppos(2)  ,  c4 ,'linewidth', 1);   % planet
    hpt = plot(ppos(1), ppos(2)  ,  c3 ,'linewidth', 1);   % planet trajectory
    scaled_a = 1.05 * a;
    axis equal                                             % ensure same scales for the two axes
    axis([-scaled_a, scaled_a, -scaled_a, scaled_a])       % set min and max axis limits
    axis manual                                            % freeze axes

    % Initializing arrays to hold values
    tsave   = zeros(1, clockmax);
    s1xsave = zeros(1, clockmax); %Star 1
    s1ysave = zeros(1, clockmax);
    s2xsave = zeros(1, clockmax); %Star 2
    s2ysave = zeros(1, clockmax);
    pxsave  = zeros(1, clockmax); %Planet
    pysave  = zeros(1, clockmax);
    dt = tmax/clockmax;
    for clock = 1:clockmax
        ind = 1;
        r = norm(X);                                  % Distance between two stars
        t = clock*dt;                                 % Updating time
        U = U - dt*G*M*X/r^3;                         % Updating velocity for star 1 and star 2
        X = X + dt*U;                                 % Updating X
    
        % Calculate individual positions and velocities from X and U
        s1pos = -m2/M*X;
        s2pos =  m1/M*X;

        % Update planet position and velocity due to star 1
        rp1          = ppos - s1pos;                    % Vector from planet to star 1
        rp1norm      = norm(rp1);                       % Distance btwn planet and star 1 
        rp2          = ppos -s2pos;                     % Vector from planet to star 2
        rp2norm      = norm(rp2);                       % Distance btwn planet and star 2
        vPS1         = G*m1*rp1/rp1norm^3;              % Force of star1 on planet
        vPS2         = G*m2*rp2/rp2norm^3;              % Force of star2 on planet
        UP           = UP - dt*vPS1 - dt*vPS2;          % Update planet's velocity due to the stars
        ppos         = ppos + dt*UP;                    % Update planet's position due to it's velocity

        % Save positions for plotting
        tsave(clock)   = t;
        s1xsave(clock) = s1pos(1);
        s1ysave(clock) = s1pos(2);
        s2xsave(clock) = s2pos(1);
        s2ysave(clock) = s2pos(2);
        pxsave(clock)  = ppos(1);
        pysave(clock)  = ppos(2);
    
        if mod(clock, 180000) == 0
            ind=ind+ind
            ht1.XData = s1xsave(1:clock);
            ht1.YData = s1ysave(1:clock);
            ht2.XData = s2xsave(1:clock);
            ht2.YData = s2ysave(1:clock);
            hpt.XData = pxsave(1:clock);
            hpt.YData = pysave(1:clock);
            hp1.XData = s1pos(1);
            hp1.YData = s1pos(2);
            hp2.XData = s2pos(1);
            hp2.YData = s2pos(2);
            hpp.XData = ppos(1);
            hpp.YData = ppos(2);
            disp(t/tmax * 100)          % Print out the percentage until completion 
            saveas(gcf, 'frame.png');
            %drawnow
        end
        hold on                         % To be able to plot multiple on same graph 
    end
en