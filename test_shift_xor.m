clear; 

a = load('vec_shift_xor_in.txt');
b = load('vec_shift_xor_out.txt');

x0 = a(:,1); y0 = a(:,2); shift_val = a(1,3); len = a(1,4);
x1 = b(:,1); y1 = b(:,2); 
x0 = x0 + circshift(y0,-mod(shift_val,len));

x0 = mod(x0,2);

err=x0-x1;

sum(err)