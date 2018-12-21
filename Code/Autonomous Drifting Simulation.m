clc; clear all; close all;
%%
p=1/2; %pick values like 1/4 loop, 1/2 loop, or 1-full loop and more
%%
% Settings for simulation time
DriftRadius = 20; % Radius of driftin meters
L=2*pi*DriftRadius; 
v=8; % speed of vehicle in m/s, this is usually applied in the kinematic bycicle model
vx=v; % longitudinal speed used in full car model.
t_loop=(L/v); %37.7 seconds is approximate time to do 1 revolution at v=5m/s
Tsim=t_loop*p; %Simulation time in seconds. the p factor gives us the number of loops.
N=20; % steps into horizon (in seconds it would be t-Horizon=N*dt)
dt=0.1; %simulation time in seconds
%Physical parameters of Car
m=2300; % in Kg
Iz=4400; % from Johns paper.
a=1.5; lf=a; %length from CoG to front axle in meters
b=1.4; lr=b; %length from CoG to back axle in meters

%%
%TIRE MODEL
Fz = (1/2)*(m*9.81)/1000; %force in K-newtons
a1 = -22.1; a2 = 1011; a3 = 1078; a4 = 1.82; a5 = 0.208; a6 = 0.000; a7 = -0.354; a8 = 0.707; C = 1.30;
D = a1*(Fz^2) + a2*Fz; BCD = a3*sin(a4*atan(a5*Fz)); B = BCD/(C*D); E = a6*(Fz^2) + a7*Fz + a8;

%%
%Declare State and Input Variables
z = sdpvar(5,N); % the oder is z=[X,Y,PSI,VY,r]'
u = sdpvar(1,N); % Delta
%declare the starting state zo. 
ztemp(:,1) = [0 -DriftRadius 0 0 0]'; % STARTING POINT

%%
for t=1:1:(Tsim/dt) 
    %Initial z value goes into constrains
    Constr = [z(:,1) == ztemp(:,t)];
    %%
    %Initial constrainst
     for i=1:N
            Constr = [Constr ];
     end
    %% ADD constraints
      
    for j=1:N-1
        alphaF = (atan( (z(4,j)+ lf*z(5,j))/vx )-u(1,j)); %WEIRD STUFF HERE
        alphaR = (atan( (z(4,j)- lr*z(5,j))/vx ));
        FyF = D*sin(C*atan(B*alphaF));
        FyR = D*sin(C*atan(B*alphaR));
                       
        Constr = [Constr,...
        %DYNAMIC EQUATIONS OF THE CAR
        z(1,j+1) == z(1,j) + (vx*cos(z(3,j))-z(4,j)*sin(z(3,j)))*dt, ... %X
        z(2,j+1) == z(2,j) + (vx*sin(z(3,j))+z(4,j)*cos(z(3,j)))*dt, ... %Y
        z(3,j+1) == z(3,j) +  z(5,j)*dt, ... %PSI
        z(4,j+1) == z(4,j) + (((FyF*cos(u(1,j))+FyR)/m)-vx*z(5,j))*dt, ... %Vy
        z(5,j+1) == z(5,j) + (((lf*FyF*cos(u(1,j))-lr*FyR)/Iz))*dt,... %r =yaw rate
        
        %% Vy constraints
        -5 <= z(4,j+1) - z(4,j) <= 5 %20m/s =72km/h
        
        %% Input constraints
        -40*pi/180 <= u(1,j+1) - u(1,j) <= 40*pi/180 
        ];
    end
    
    %% CREATE REFRENCE
    z_ref = getref(ztemp(1,t),ztemp(2,t),dt,N,DriftRadius,vx);
        
    %%
    % COST FUNCTION
    Cost = 0;
    for k=1:N
        Cost = Cost + norm( z(1:2,k)-z_ref(1:2,k) )^2; % we are only using the X_ref and Y_ref
    end
    %%
    %Solve Optimization Problem
    options = sdpsettings('verbose', 1,'solver','ipopt');
    solvesdp(Constr,Cost,options)
    %display a timer for ease of reading
    disp(['t= ',num2str(t*dt),' of ',num2str(Tsim),' seconds in increments of ',num2str(dt),' seconds']); %Allows you for visual timer
    %%
    %use the optmized U of the solver and save them.
    u_star(:,t) = double(u(:,1)); %Use only the first value of the optimized solution
    zopen{1,t} = double(z);
        
    %Update the Alfas and Forces.
    alphaF(t) = double( atan( (ztemp(4,t)+ lf*ztemp(5,t))/vx )-u_star(1,t) );
    alphaR(t) = double( atan( (ztemp(4,t)- lr*ztemp(5,t))/vx ) );
    FyF(t) = D*sin(C*atan(B*alphaF(t))); %trick to include the mass of the car FIX FIX FIX
    FyR(t) = D*sin(C*atan(B*alphaR(t))); % trick to include the mass of th car FIX FIX FIX
    % the first value of the optimized U + the z to predict the next z.
    ztemp(1,t+1) = ztemp(1,t) + (vx*cos(ztemp(3,t))-ztemp(4,t)*sin(ztemp(3,t)))*dt; %X
    ztemp(2,t+1) = ztemp(2,t) + (vx*sin(ztemp(3,t))+ztemp(4,t)*cos(ztemp(3,t)))*dt; %Y
    ztemp(3,t+1) = ztemp(3,t) + ztemp(5,t)*dt; % PSI
    ztemp(4,t+1) = ztemp(4,t) + (((FyF(t)*cos(u_star(1,t))+FyR(t))/m)-vx*ztemp(5,t))*dt; %Vy
    ztemp(5,t+1) = ztemp(5,t) + (((lf*FyF(t)*cos(u_star(1,t))-lr*FyR(t))/Iz)*dt); %r
    
    %%
    %PLOT OF X VS Y, SIMULATION MPC VS REFERENCE
    plot(ztemp(1,:),ztemp(2,:),'b--o','LineWidth',2); hold on; %Plot open loop
    circle(0,0,DriftRadius) % Plot Circl trajectory
    plot(z_ref(1,:),z_ref(2,:),'r--o'); hold off;
    title('X vs Y MPC SIMULATION');
    %%axis([-10,60,-52,10]);
    
    %%
end

%plot Trajectory
plot(ztemp(1,:),ztemp(2,:),'b--o','LineWidth',2);
axis([-60,60,-60,60])

%Drift Animation
%CarTrajectoryPlot([1:t],ztemp',u_star',-DriftRadius-20,DriftRadius+20,-DriftRadius-20,DriftRadius+20,[])

%Vy
figure();plot(ztemp(4,:))

%PLOT INPUT
figure()
plot(dt*[1:1:length(u_star)],(180/pi)*u_star(1,:),'b--o'); hold on; %DELTA=[Degrees]
title('input DELTA'); legend('DELTA');xlabel('seconds'); ylabel('degrees');