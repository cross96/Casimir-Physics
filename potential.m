polarizabilities = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM1_ESM.xlsx");
polarizabilities1 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM2_ESM.xlsx");
polarizabilities2 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM3_ESM.xlsx");
polarizabilities3 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM4_ESM.xlsx");
polarizabilities4 = readtable("C:\Users\cross\Desktop\polarizability_tables\10053_2013_529_MOESM5_ESM.xlsx");

kb = 1.3807e-23;
h_bar = 6.62607015e-34;
c = 299792458;
T = 10;

points = [-0.0185500000000000,
-0.00927500000000000,
0,
0.00927500000000000,
0.0185500000000000,
0.0278250000000000,
0.0371000000000000,
0.0463750000000000,
0.0556500000000000,
0.0649250000000000,
0.0742000000000000]

y' = wick_rotate(FEL(G));

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);


t=-3.5*st:dt:3.5*st;
m = morlet_m(f,t,width);
y = conv(s,m);

C = [0+i 1+i 1+100i 0+100i];
q = integral(fun Fs*sf/st, 1,1, C, C);
%TODO add waypoints
%TODO Eignevalues
%check = 0;
%while check==0
%   clear N ; 
%   for N = 1:MaxMod ;
%       % Data
%       p = modord(N,Ni) ; % p represents the order of the linear differential equations 
%       % Overdetermined matrix definition
%%       % Overdetermined equation A*x=B formulation
 %      B = -G(:,1:Ni) ;
 %      A = G(:,Ni+1:(p+1)*Ni) ;
 %      % Resolving (using the pseudoinverse function)
 %      x = pinv(A)*B ; 
 %      % Erreur calculation 
 %      %		- Inverse of conditionnning number
 %      InvCond(N) = 1/cond(A) ;
 %%      [U,S,V] = svd(A,0) ; 
  %     SingVals = diag(S) ;
  %     err(N) = 1/(max(SingVals)/min(SingVals)) ;
  %     %		- of least squares
  %     epsilon(N) = norm(B-(A*x));
  %     % Eigenvalue problem resolving
  %     clear L z ;
  %     [L,z] = PbValPp(x,Ni,p) ;
  %     % Modal parameters extraction
  %     clear lambda wd delta wn xi i ;
  %     lambda = log(z)./deltaT ;
  %     wd = imag(lambda) ;
   %    delta = real(lambda) ;
   %    wn = sqrt(wd.^2+delta.^2) ;
   %    xi = -(delta./wn) ;
   %    fn = wn/(2*pi) ;
       
function wick_rotate = e_c(omega, r)
e_c = omega^2*i/r * lookup_table(omega);
end
function lookup_table = alpha(omega)
lookup_table = interp1(omega, polarizabilities, nearest);
end

function Ej_Ek = FEL(G)
Ej_Ek = zeros(2,100);
for w_r = 1+250e12:100+250e12
    for w_i = 1+250e12:100+250e12
        Ej_Ek(w_r-250e12,w_i-250e12) = -6.62607015e-34/pi .* imag((w_r +w_i).^2 .* G(w_r-250e12)).*coth(6.62607015e-34 .* (w_r+w_i) ./ 2*1.3807e-23);
    end
end
end
