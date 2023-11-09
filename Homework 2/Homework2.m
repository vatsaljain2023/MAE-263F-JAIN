
% DOF vector:
% q = [x1; y1; z1; theta1; x2; y2; z2; theta2; ...
% theta_{N-1}; xN; yN; zN];
% nv = number of vertices/nodes
% ne = nv - 1: number of edges
% ndof = size(q) = 4*nv-1 = 3*nv + ne
% Problem: q(t) = ? Given Initial condition (q(0))
% and physical parameters
%
% function tangent = computeTangent(q)
%
% function d = parallel_transport(u, t1, t2)
%
% Time parallel reference frame:
% function [a1, a2] = computeTimeParallel(a1_old, q0, q)
%
% Material frame: function [m1,m2] = ...
% computeMaterialDirectors(a1, a2, theta)
%
% function refTwist = computeRefTwist(a1, tangent, refTwist)
%
close all;
clear all;
clc

%% global variables
global Fg M dt
global kappaBar EI GJ voronoiLength
global EA refLen

%% Inputs

% (1) Define the number of DOF related quantities
nv = 50; % number of nodes or vertices
ne = nv - 1; % number of edges
ndof = 3*nv + ne; % number od DOG

% (2) Define time step related quantities
dt = 0.01; % time step size in second
totalTime = 5; % simulation time
Nsteps = round(totalTime / dt);

% (3) Physical parameters
% (3a) Geometry
RodLength = 0.2; % meter
natR = 0.02; % meter
r0 = 0.001; % meter

% (3b) Material parameters
Y = 10e6; % Young's modulus
nu = 0.5; % Poisson's ratio
G = Y / (2 * (1+nu)); % Shear modulus
rho = 1000; % density (kg/m^3)

% (3c) Geometry
g = [0;0;-9.81];

%% Stifness variables
EI = Y * pi * r0^4 / 4; % Bending stiffness
GJ = G * pi * r0^4/2; % Shearing stiffness
EA = Y * pi * r0^2; % Stretching stiffness

%% Tolerance
tol = EI / RodLength^2 * 1e-6;

%% Mass Matrix
totalM = pi*r0^2*RodLength*rho; % kg
dm = totalM / ne; % mass per edge
massVector = zeros(ndof, 1);

for c = 1:nv
    ind = [4*c-3; 4*c-2; 4*c-1]; % c-th node
    if c==1 || c==nv
        massVector(ind) = dm/2;
    else
    massVector(ind) = dm;
    end
end

for c=1:ne
    ind = 4*c;
    massVector(ind) = 1/2 * dm * r0^2;
end

M = diag(massVector);

%% Initial DOF vector
nodes = zeros(nv, 3);
dTheta = (RodLength / natR) * (1/ne);

for c=1:nv
    nodes(c,1) = natR * cos ( (c-1) * dTheta );
    nodes(c,2) = natR * sin ( (c-1) * dTheta );
    nodes(c,3) = 0;
end

q0 = zeros(ndof, 1); % Initial condition

for c=1:nv
    ind = [4*c-3; 4*c-2; 4*c-1]; % c-th node
    q0(ind) = nodes(c,:);
end

% All theta angles are zero, i.e., q0(4:4:end) = 0 but it's redundant
u = zeros(ndof, 1); % Velocity

%% Reference length for each edge (used in stretching force)
refLen = zeros(ne, 1);
for c=1:ne % loop over each edge
    dx = nodes(c+1,:) - nodes(c,:); % dx = x_{c+1} - x_c
    refLen(c) = norm(dx);
end

%% Voronoi length (length associated with each node; used for bending and twisting)
voronoiLength = zeros(nv, 1);

for c=1:nv
    if c==1
        voronoiLength(c) = 1/2 * refLen(c);
    elseif c==nv
        voronoiLength(c) = 1/2 * refLen(c-1);
    else
        voronoiLength(c) = 1/2 * refLen(c-1) + 1/2 * refLen(c);
    end
end

%% Reference frame (At t=0, initialize using space parallel transport)
a1 = zeros(ne, 3); % First reference director for all the edges
a2 = zeros(ne, 3); % Second reference director for all the edges
tangent = computeTangent( q0 ); % Tangent for all the edges; size is (ne,3)

% Compute a1 for the first edge
t0 = tangent(1, :); % tangent on the first edge
t1 = [0;0;-1]; % arbitrary
a1Tmp = cross(t0, t1); % This vector is perp. to t0

if abs(a1Tmp) < 1e-6
    t1 = [0;1;0]; % arbitrary
    a1Tmp = cross(t0, t1); % This vector is perp. to t0
end

a1(1,:) = a1Tmp / norm(a1Tmp);
a2(1,:) = cross(tangent(1,:), a1(1,:));

% Done with the first edge
for c=2:ne
    t0 = tangent(c-1,:); % previous tangent on c-1-th edge
    t1 = tangent(c,:); % current tangent on c-th edge
    a1_0 = a1(c-1,:); % previous a1 director on c-1-th edge
    a1_1 = parallel_transport( a1_0, t0, t1);
    a1(c,:) = a1_1 / norm( a1_1 );
     a2(c,:) = cross(t1, a1(c,:));
end

%% Material frame
theta = q0(4:4:end); % Vector of size ne
[m1, m2] = computeMaterialDirectors( a1, a2, theta);

%% Reference twist
refTwist = zeros(nv, 1);

%% Natural curvature
kappaBar = getkappa( q0, m1, m2 ); % Natural curvature at each node

%% Gravity
Fg = zeros(ndof, 1); % Gravity vector

for c=1:nv % loop over the nodes
    ind = [4*c-3; 4*c-2; 4*c-1]; % c-th node
    Fg(ind) = massVector(ind) .* g; % vector of size 3
end

%% Fixed and free DOFs
fixedIndex = 1:7; % Clamped BC
freeIndex = 8:ndof;

%% Time stepping scheme
ctime = 0; % Current time
endZ = zeros(Nsteps, 1); % z-coordinate of the last node with time

for timeStep = 1:Nsteps
    fprintf('Current time = %f\n', ctime);
    [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, ...
    tol, refTwist);
    ctime = ctime + dt;

    % Update q0 (old position)
    q0 = q;

    % Store endZ
    endZ(timeStep) = q(end);

    % Plot
    if mod(timeStep,100) == 0
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors( a1, a2, theta);
        plotrod(q, a1, a2, m1, m2, ctime);
    end
end

% Visualization
figure(2);
timeArray = (1:1:Nsteps) * dt;
plot( timeArray, endZ, 'ro-');
xlabel('Time, t [sec]');
ylabel('z-coordinate of last node, \delta_z [m]');


