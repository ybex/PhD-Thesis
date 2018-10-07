function [Sx, Sy, Sz] = make_spin_matrices_l(S)

Sx = zeros(2*S+1);
Sy = zeros(2*S+1);
Sz = zeros(2*S+1);
Splus = zeros(2*S+1);
Sminus = zeros(2*S+1);

for i = 1:2*S 
    m = -S+(i-1); 
    Splus(i+1,i) = sqrt(S*(S+1)-m*(m+1)); 
end
for i = 2:2*S+1 
    m = -S+(i-1); 
    Sminus(i-1,i) = sqrt(S*(S+1)-m*(m-1)); 
end
Sx = (Splus + Sminus)/2;
Sy = -1i * (Splus - Sminus)/2;
for i = 1:2*S+1 
    Sz(i,i) = -S+(i-1); 
end