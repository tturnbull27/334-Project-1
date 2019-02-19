function [g_new] = Turnbull_1002243921_g_values(td,n)

g_new(1,1) = td(1,2)/td(1,1);

for i = 2 : n-1
    
    g_new (i,1) = td(i,i+1)/(td(i,i)-(td(i-1,i)*g_new(i-1,1)));

end
