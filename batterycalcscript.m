Checkpoints=10;
Distance=1000;
alpha=20;
% alpha in degrees
mass=17;
% mass in kg
Cd=1.3;
g=9.8066;
% g in m/(s^2)
wind_velocity=[0 5];
%  wind velocity in m/s
area=1.55;
% area in m^2
T=mass*g/cosd(alpha)/6;
% Maximum achievable Thrust
fprintf('Maximum achievable Thrust per motor = %.4f N\n',T);
T_sin_alpha_by_m=g*tand(alpha);
constant_for_drag=0.5*1.225*area*Cd;
% Density of air=1.225 kg/m^3
K=constant_for_drag*sind(alpha)*(1/mass);

Battery=[0 0];
Total_time=[0 0];

%I=[3.08,3.45,3.87,4.39,4.89,5.43,6.06,6.71,7.35,8.01,8.68,9.37,10.22,10.99,11.73,12.62,14.75,17.13,22.59,27.75];
%thrust=[889,968,1058,1172,1271,1382,1489,1592,1712,1829,1944,2078,2163,2310,2424,2535,2852,3068,3710,4284]*g/1000;
I=[3.45,4.39,5.43,6.71,8.01,9.37,10.99,12.62,14.75,17.13,22.59,27.75];
thrust=[968,1172,1382,1592,1829,2078,2310,2535,2852,3068,3710,4284]*g/1000;
[const,e,mu]=polyfit(thrust,I,8);
Current=polyval(const,T,e,mu);
fprintf('Current required = %.4f amperes\n',Current);

for n1=1:2
Drone_max_velocity=((mass*g/(cosd(alpha)*constant_for_drag))^0.5)-wind_velocity(n1);    

fprintf('--------------------------------------------------------------------\n');
fprintf('Wind Velocity = %d m/s\n',wind_velocity(n1));

syms x
velocity_start=0;
velocity_end=0;
diff_time_Acc=(1/(T_sin_alpha_by_m-(K*((wind_velocity(n1)+x)^2))));
diff_time_Dec=-(1/(T_sin_alpha_by_m+(K*((wind_velocity(n1)+x)^2))));
% sind(alpha) term is for component of area
time_Acc=int(diff_time_Acc,velocity_start,Drone_max_velocity);
time_Dec=int(diff_time_Dec,Drone_max_velocity,velocity_end);

syms s(t)
Ds = diff(s);
Disp_function_Acc = diff(s,t,2)+(K*(diff(s,t,1))^2)+(2*K*wind_velocity(n1)*diff(s,t,1)) == T_sin_alpha_by_m-((wind_velocity(n1)^2)*K);
Disp_function_Dec = -diff(s,t,2)-(K*(diff(s,t,1))^2)-(2*K*wind_velocity(n1)*diff(s,t,1)) == T_sin_alpha_by_m+((wind_velocity(n1)^2)*K);
cond1 = s(0) == 0;
cond2 = Ds(0) == 0;
cond3 = Ds(0) == Drone_max_velocity;
conds1 = [cond1 cond2];
conds2 = [cond1 cond3];
Disp_Acc(t) = dsolve(Disp_function_Acc,conds1);
Disp_Dec(t) = dsolve(Disp_function_Dec,conds2);

t=time_Acc;
d_Acc=double(Disp_Acc(t));
t=time_Dec;
d_Dec=double(Disp_Dec(t));
D=d_Acc+d_Dec;

if Distance > D
    time_const_velocity=(Distance-D)/Drone_max_velocity;
    Total_time(n1) = time_const_velocity + time_Acc + time_Dec;
    fprintf('max velocity drone can achieve = %.4f m/s\n',Drone_max_velocity);
    fprintf('time for acceleration = %.4f seconds\n',time_Acc);
    fprintf('time for deceleration = %.4f seconds\n',time_Dec);
    fprintf('distance covered for acceleration = %.4f meters\n',d_Acc);
    fprintf('distance covered for constant velocity = %.4f meters\n',Distance-D);
    fprintf('distance covered for deceleration = %.4f meters\n',d_Dec);
else
    syms v
    Diff_eqn_Disp_function_Acc = (v/(T_sin_alpha_by_m-(K*((wind_velocity(n1)+v)^2))));
    Diff_eqn_Disp_function_Dec = -(v/(T_sin_alpha_by_m+(K*((wind_velocity(n1)+v)^2))));
    Disp_function_Acc(v) = int(Diff_eqn_Disp_function_Acc);
    Disp_function_Dec(v) = int(Diff_eqn_Disp_function_Dec);
    Disp_function_Acc(v) = Disp_function_Acc(v) - Disp_function_Acc(0);
    Disp_function_Dec(v) = Disp_function_Dec(0) - Disp_function_Dec(v);

    Eq1 =  Disp_function_Acc(v) + Disp_function_Dec(v) == Distance;
    velocity=vpasolve(Eq1,[0 Inf]);
    Disp_Acc = int(Diff_eqn_Disp_function_Acc,[0 velocity]);
    Disp_Dec = int(Diff_eqn_Disp_function_Dec,[velocity 0]);
   
    syms x
    time_Acc=int(diff_time_Acc,0,velocity);
    time_Dec=int(diff_time_Dec,velocity,0);
    Total_time(n1) = time_Acc + time_Dec;
   
    fprintf('max velocity drone can achieve = %.4f m/s\n',Drone_max_velocity);
    fprintf('time for acceleration = %.4f seconds\n',time_Acc);
    fprintf('time for deceleration = %.4f seconds\n',time_Dec);
    fprintf('distance covered for acceleration = %.4f meters\n',Disp_Acc);
    fprintf('distance covered for deceleration = %.4f meters\n',Disp_Dec);
end

Battery(n1)=(Current*1000)*(Total_time(n1)/3600)*6*(Checkpoints/2)/0.8;
% Current in amperes and Total_time in seconds, no. of motors = 6
% factor of 1/0.8 considering 80% efficiency
end
fprintf('--------------------------------------------------------------------\n');

Battery=sum(Battery);
Total_time=sum(Total_time);
fprintf('Total time for %d checkpoints = %.4f minutes\n',Checkpoints,(Checkpoints/2)*Total_time/60);
fprintf('Total battery capacity required for %d checkpoints = %.4f mAh\n',Checkpoints,Battery);
