function [r_new] = Turnbull_1002243921_r_values(td,r,g_new,n)

r_new(1,1)= r(1,1)/td(1,1);

for i=2: n
    
    r_new(i,1) = (r(i,1)-td(i,i-1)*r_new(i-1,1))/(td(i,i)-td(i,i-1)*g_new(i-1,1));

end 

    