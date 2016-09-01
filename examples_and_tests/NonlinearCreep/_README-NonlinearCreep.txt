
This file explains the contents of the CHILD test repository folder
called "NonlinearCreep".

The folder contains a series of test input files and selected output
files that demonstrate the performance of the nonlinear creep algorithm
as implemented in CHILD as of August, 2007.

There are four test cases. All parameters are identical in each of these
except the uplift rate, which varies from 10^-5 to 10^-3 m/yr:

Test:                     Uplift rate (m/yr):
-----                     -------------------

nldifftestbasin1s         0.001
nldifftestbasin2s         0.00001
nldifftestbasin3s         0.0001
nldifftestbasin4s         0.0003

Relief increases systematically with uplift rate, of course, and as it
does, the importance of the nonlinear portion of the qs vs. S curve
grows. I don't know that these runs are strictly at steady state by
t=3e6 (1e6 in the case of 1s) but they appear to be close based on
slope-area plots. Note that run 1s is a lot slower because most of the
terrain is well within the nonlinear realm and thus requires very small
time steps.

Note that a separate series of tests, with only creep active, were also
run and compared with analytical solutions. I have not stored these 
files. However, numerical solutions to linear and nonlinear creep can 
be compared easily by setting up a domain of, say, 10m by 100m with
one open-side boundary. The equilibrium solution for linear creep is:

     U     2    2
z = ---- (L  - x  )
    2 Kd

where L is the length of the domain (e.g., 100m) and x = L - y, where y 
is CHILD's y coordinate (=0 at the boundary).

For nonlinear creep, there is an analytical solution for equilibrium slope
angle that has two roots. Experimentally, it seems that the first root
is the one that matters, and it looks like:

                                         2      2 1/2
                               Sc (a - (a  + 4 b )   )
                     S = - 1/2 -----------------------
                                          b

where Sc is critical slope, a = Kd Sc, and b = U x. Beware however that
CHILD's *.slp file contains the slopes of the triangle edges, which are
generally not parallel to the y-axis. So, the curve above should be an
upper envelope for slope (though numerical error will undoubtedly leave
some points slightly above). To do it properly, one could retrieve the
y-component of slope for each edge, but I haven't tried to do this.

The Matlab script below uses the same algorithm in CHILD (as of Aug 07
anyway) to solve the nonlinear creep equation in 1d. (The equation is that
of Howard, 1994).


%NLDIFF1D: simple model of a hillslope profile with nonlinear creep
%transport.
%
% Left boundary condition is no flux, representing a drainage divide.
% Right boundary condition is specified elevation, representing baselevel.
%
% The time step size is based on a sort of linearization. The equation is
%            S
% qs = kd ----------, combined with mass continuity
%         (1-|S|/Sc)
%
% With this, I estimate a Courant condition:
% dt <= dx^2 / 2 (kd / (1-(|S|/Sc)^2))
%
% But I find instability as |S|->Sc. This instability seems to go away if I
% multiply by an extra sqrt(1-(|S|/Sc)^2).
%
% Also, note that in order to prevent the equation blowing up or going
% seriously askew with |S|/Sc >= 1, I force an upper bound (beta), which
% can effect solutions when they are very steep with a high flux.
%
%GT, Aug 2007
%

% Parameters and vectors
L = 100;   % hillslope length, meters
nn = 100;  % number of nodes
dx = L/nn; % node spacing
xc = (dx/2):dx:L;   % cell centers, including right-hand boundary cell
kd = 1e-3;  % creep coefficient, m2/yr
u = 1e-3;   % baselevel fall rate, m/yr
remtm = 1e5;  % run time (and remaining time), yr
z = zeros(1,nn)+u*remtm;   % initial elevations
Sc = 0.8;   % threshold slope gradient
beta = 0.999;  % effective maximum S/Sc
xs = 0:dx:(L-dx);  % positions of interfaces between nodes, incl left edge (0)

% Plotting stuff
plotinterval = 0.1*remtm;   % plot interval, years
nextplot = plotinterval;
eltm = 0;   % elapsed time in run, yr
figure(1)
clf
subplot(2,1,1)
plot(xc,z)
axis equal
axis([0 L 0 1.2*max(z)])
hold on
subplot(2,1,2)
plot(xs,zeros(length(xs)));
ylabel('Slope')
hold on
drawnow

% Main loop
while remtm>1e-6
   
    S = [0 -diff(z)/dx];  % calculate slopes between nodes (0 is left bdy)
    phi = min( abs(S)/Sc, beta );  % ratio of slope to critical slope
    fluxcoef = kd./(1-phi.^2);    % "flux coefficient", ie what's multiplied by S
    dts = 0.1*dx^2.*sqrt(1-phi.^2)./fluxcoef;  % maximum time step for each node
    dt = min( dts );   % smallest of these
    dt = min( dt, remtm );  % don't go past the end of the run
    if dt<=0.0   % this should never happen
        error('Zero or negative dt');
    end
    qs = fluxcoef.*S;   % transport rates (m2/yr)
    dqsdx = diff(qs)/dx;  % gradient in transport rates
    dz = -dqsdx*dt;       % change in elevation this iteration
    z(1:nn-1) = z(1:nn-1) + dz;  % update elevations
    z(nn) = z(nn) - u*dt;   % lower baselevel
    remtm = remtm - dt;
    eltm = eltm + dt;
    
    if eltm>=nextplot
       subplot(2,1,1)
       plot(xc,z)
       subplot(2,1,2)
       plot(xs,S)
       nextplot = nextplot+plotinterval;
       drawnow
    end
    
end


% plot analytical solution for steady state, linear case
sa = xs*(u/kd);
subplot(2,1,2)
hold on
plot(xs,sa,'ro')

% plot analytical solution for steady state, nonlinear
% Note: the quadratic expressions are the two roots for V = S/Sc, so we
% multiply each set of roots by Sc to recover S.
a = kd*Sc;
b = u*xs;
snl1= -1/2*(a-(a^2+4*b.^2).^(1/2))./b;
snl2 = -1/2*(a+(a^2+4*b.^2).^(1/2))./b;
snl1 = snl1*Sc;
snl2 = snl2*Sc;
plot(xs,snl1,'go')
plot(xs,snl2,'ko')

% Clean up
hold off
subplot(2,1,1)
hold off

