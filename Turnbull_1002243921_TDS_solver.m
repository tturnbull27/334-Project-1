
function [x] = Turnbull_1002243921_TDS_solver(D,r)

%this function solves the TDS matrix and returns the unknown: x
%td is the tridiagonal matrix, and r is the matrix of constants

%find the size of td and put it in to n for for loops
size_td = size(D);
n = size_td(1,1);

%call the function to solve for the gi values
[g_new] = Turnbull_1002243921_g_values(D,n);

% call the function to modify the constants matrix
[r_new] = Turnbull_1002243921_r_values(D,r,g_new,n);

%calculating the values using backsub
x(n,1) = r_new(n,1);

for i = 2:n
    
    j = n-i+1;
    
    x(j,1) = r_new(j,1)-g_new(j,1)*x(j+1,1);
    
end







