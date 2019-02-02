polarizabilities1 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM1_ESM.xlsx");
%polarizabilities2 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM2_ESM.xlsx");
%polarizabilities3 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM3_ESM.xlsx");
%polarizabilities4 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM4_ESM.xlsx");
%polarizabilities5 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM5_ESM.xlsx");

kb = 1.3807e-23;
h_bar = 6.62607015e-34;
c = 299792458;
T = 10;
e = 1;
m = 1;
w = 400;
L = 10e-9;
e0 = 8.85418782e-12;

 E_mag = readtable("C:\Users\cross\Desktop\E_mag.txt");
 E_y_Re = readtable("C:\Users\cross\Desktop\E_y_Re.txt");
 E_y_Im = readtable("C:\Users\cross\Desktop\E_y_Re.txt");

 greens_fun_x = @(k) k*(k^2-w^2/c^2)^0.5*(((e*(k^2-w^2/c^2)^.5...
     -(k^2-e*m*w^2/c^2)^.5))/(e*(k^2-w^2/c^2)^.5...
     +(k^2+e*m*w^2/c^2)^.5) +w^2/(c^2*(k^2-w^2/c^2)^.5)*((k^2-w^2/c^2)^...
     0.5*(((m*(k^2-w^2/c^2)^.5...
     -(k^2-e*m*w^2/c^2)^.5))/(m*(k^2-w^2/c^2)^.5...
     +(k^2+e*m*w^2/c^2)^.5))))*exp(-2*k*L);
  
 greens_fun_y = greens_fun_x;
  
 greens_fun_z = @(k) k *2*k^2/(k^2-w^2/c^2)*(e*(k^2-w^2/c^2)^.5-...
     (k^2-e*u*w^2/c^2)^.5)/(e*(k^2-w^2/c^2)^.5+(k^2-e*u*w^2/c^2)^.5);
 
 