
function[T, x] = Turnbull_1002243921_project_1(properties,dimensions,n)

hold on;
w = 1;
L = dimensions(1,1);
b = dimensions(2,1);
q_b = properties(1,1);
T_inf = properties(2,1);
h_max = properties(3,1);
h_min = properties(4,1);
k = properties(5,1);

%Useful parameters
theta = atan(b/(2*L));
delta_x = L/(n-1);
A_b = b*w;

%Initialize matrices
D = zeros(n,n);
r = zeros([n 1]);
x = zeros([n 1]);


% Boundary condition at base of fin
i = 1;
a_right= 2*w*(L-(i-0.5)*delta_x)*tan(theta);
a_left= 2*w*(L-(i-1.5)*delta_x)*tan(theta);
a_c= delta_x/cos(theta);
h_x = (-(h_max - h_min)/L)*delta_x*i + ((h_max - h_min)/L)*delta_x + h_max;

D(1,1) = (-k * a_right) / delta_x  -  a_c * h_x ;
D(1,2) = (k * a_right) / delta_x ;
r(1,1) = (-T_inf) * a_c * h_x  -  q_b * A_b ;

%Boundary condition at node n
i = n;
a_right= 2*w*(L-(i-0.5)*delta_x)*tan(theta);
a_left= 2*w*(L-(i-1.5)*delta_x)*tan(theta);
a_c= delta_x/cos(theta);
h_x = (-(h_max - h_min)/L)*delta_x*i + ((h_max - h_min)/L)*delta_x + h_max;

D(n,n)   = (-k * a_left) / delta_x  -  2*h_x*a_c ;
D(n,n-1) = ( k * a_left) / delta_x ;
r(n,1)   = (-2)*h_x*a_c*T_inf ;
x(n,1)   = L ;

%Filling in the rest of the tridiagonal matrix, r matrix, and x matrix.  
for i=2:(n-1)
    
    a_right= 2*w*(L-(i-0.5)*delta_x)*tan(theta);
	a_left= 2*w*(L-(i-1.5)*delta_x)*tan(theta);
    a_c= delta_x/cos(theta);
    h_x = (-(h_max - h_min)/L)*delta_x*i + ((h_max - h_min)/L)*delta_x + h_max;
    
    D(i,i-1) = ( k*a_left)/delta_x ;
    D(i,i)   = (-k*a_left)/delta_x - (k*a_right)/delta_x - 2*h_x*a_c;
    D(i,i+1) = ( k*a_right)/delta_x ;
    r(i,1)   = (-2) * h_x * a_c * T_inf ;
    x(i,1)   = delta_x*(i-1);
    
end

[T] = Turnbull_1002243921_TDS_solver(D,r);  %Array of temperature values for each node

plot(x,T);
title("Temperature Distribution for Various Control Volumes");
xlabel("x");
ylabel("Temperature");
legend('Variable h', 'Constant h')

end