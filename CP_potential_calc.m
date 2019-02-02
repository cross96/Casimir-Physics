%Import data from Lumerical as needed)

polarizabilities1 = xlsread("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM1_ESM.xlsx");
%polarizabilities2 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM2_ESM.xlsx");
%polarizabilities3 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM3_ESM.xlsx");
%polarizabilities4 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM4_ESM.xlsx");
%polarizabilities5 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM5_ESM.xlsx");
% E_mag = readtable("C:\Users\cross\Desktop\E_mag.txt");
%E_y_Re = readtable("C:\Users\cross\Desktop\E_y_Re.txt");
% E_y_Im = readtable("C:\Users\cross\Desktop\E_y_Re.txt");
%eps = readtable("C:\Users\cross\Desktop\esp.txt");

% Wick Rotations
for j = 1:(polarizabilities1(:,1))
   ec =  j^2* e / 1.0i;
end
kb = 1.3807e-23;
h_bar = 6.62607015e-34;
Pi = 3.14159265358979;
c = 299792458;
T = 10;
m = 1;
w = 400;
e0 = 8.85418782e-12;
e = 2*e0;

for L = [1.436e-07,1.343e-07,1.25e-07,1.157e-07,1.1314e-07,9.72e-08,8.79e-08,-3.39e-07,6.93e-08,6.01e-08,5.08e-08]
     greens_fun_x = @(k)@(w) k*(k^2-w^2/c^2)^0.5*(((e*(k^2-w^2/c^2)^.5...
         -(k^2-e*m*w^2/c^2)^.5))/(e*(k^2-w^2/c^2)^.5...
         +(k^2+e*m*w^2/c^2)^.5) +w^2/(c^2*(k^2-w^2/c^2)^.5)*((k^2-w^2/c^2)^...
         0.5*(((m*(k^2-w^2/c^2)^.5...
         -(k^2-e*m*w^2/c^2)^.5))/(m*(k^2-w^2/c^2)^.5...
         +(k^2+e*m*w^2/c^2)^.5))))*exp(-2*k*L);

     greens_fun_y = greens_fun_x;

     greens_fun_z = @(k)@(w) k *2*k^2/(k^2-w^2/c^2)*(e*(k^2-w^2/c^2)^.5-...
         (k^2-e*u*w^2/c^2)^.5)/(e*(k^2-w^2/c^2)^.5+(k^2-e*u*w^2/c^2)^.5)*...
         exp(-2*k*L);
    greens_tensor_x = 1 / (8 * Pi * e0) * integral(greens_fun_x,0, Inf);
 
    greens_tensor_y = greens_tensor_x;
 
    greens_tensor_z = 1 / (8 * Pi * e0) * integral (greens_fun_y,0,Inf);
    
    U = -kb * symsum(greens_tensor_x+greens_tensor_y+greens_tensor_z,400,1600);
end