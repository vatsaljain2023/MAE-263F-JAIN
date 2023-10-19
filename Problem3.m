%% Related to Assignment 1
% Number of nodes
N = 50;
ndof = N * 2; % number of degrees of freedom
dt = 0.01; % second - Time step size
totalTime = 1; % second - total simulation time
RodLength = 1; % meter
deltaL = RodLength / (N-1);

% Density
% kg/m^3
rho = 2700; %aluminum
Y = 70e9; % Young's modulus (Y instead of E for clarity)

% Radii of spheres
Ro = 0.013; %outer radius
Ri = 0.011; %inner radius

% Utility parameter
EI = Y * pi/4 * (Ro^4 - Ri^4); % Nm^2 - bending stiffness
EA = Y * pi * (Ro^2-Ri^2);

% Geometry - initial configuration
nodes = zeros(N,2);

for c=1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);

for k=1:N
    M(2*k-1, 2*k-1) = pi*(Ro^2 - Ri^2)*rho*RodLength / (N-1); % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Force vector, P
P = zeros(ndof, 1);
Pnode = 38; % 0.755 away from left edge
P(Pnode*2)=-20000; % y component of Pnode

% Initial DOF
q0 = zeros(ndof, 1);

for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end

u0 = zeros(ndof, 1); % old velocity (initial velocity)

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected

% Time marching scheme
Nsteps = round(totalTime/dt);

%storage for max vertical displacement
ymax = zeros(Nsteps,1);

for c = 2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt);
    
    % Fixed and free DOFs
    fixedDOF = [1;2; ndof];
    freeDOF = 3:ndof-1;
    boundaryConditionVector = [0;0;0]; % [0;0];

    % Guess
    q = q0; % New DOFs are initialized to be equal to old DOFs
    q(fixedDOF) = boundaryConditionVector; % IMPOSE
   
    % Newton Raphson
    err = 10 * tol;

    while err > tol
        q_free = q(freeDOF);
        
        f = M / dt * ( (q-q0)/dt - u0 );

        J = M / dt^2;

        % ELASTIC FORCES

        % Linear spring
       
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end
        
        % Bending spring
        
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            l_k = deltaL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end

        % Force
        f = f - P;

        % At this point, we have f and J
        f_free = f(freeDOF);
        J_free = J(freeDOF, freeDOF);

        % Update
        dq_free = J_free \ f_free;
        q_free = q_free - dq_free;

        err = sum ( abs(f_free) );
        q(freeDOF) = q_free;
        q(fixedDOF) = boundaryConditionVector;
    end

    % New velocity
    u = (q - q0) / dt;

    % Plot
    figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis equal
    xlabel('x [meter]');
    ylabel('y [meter]');
    drawnow
    
    %store y coordinates for this time
    ycoord=zeros(N,1);
    for i = 1:1:N
        ycoord(i) = q(2*i);
    end

    %store max displacement
    ymax(c) = max(abs(ycoord));

    % Update (new becomes old)
    q0 = q;
    u0 = u;
end

% Plot max displacement
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, -1*ymax, 'k-');
xlabel('Time, t [sec]');
ylabel('Max Displacement (vertical) [m]');
title('Max Vertical Displacement of Beam vs. Time')



