clc; clear; close
Nplies = input('Write number of layer'); %In our project all the laminates have 4 layer
h_ply  = input('Write layer thickness value with in mm =');  % for thickness of layer   
h  = Nplies * h_ply ; %For laminate thickness
for i = 1:Nplies;
     thetadt(i) = input(' Write the layer orientation from first to last =');
end
thetadb = fliplr(thetadt);

% Entering with MPa because MPa = N/mm
E_1 = input(' Write Modulus in fibre direction E_1 value with MPa =');%Modulus in fibre direction
E_2 = input(' Write Modulus transverse to fibre E_2 value with MPa ='); %Modulus transverse to fibre
G_12 = input(' Write Shear moduli G_12 value with MPa =');%Shear moduli
G_13 = input(' Write Shear moduli G_13 value with MPa =');
G_23 = input(' Write Shear moduli G_23 value with MPa =');
v_12 = input(' Write  Poisson’s ratio v_12 value ='); % Poisson’s ratio:
v_21 = v_12 * (E_2/E_1) ;
v_23 = input(' Write  Poisson’s ratio v_23 value =');

%For different load cases I used input function to use in-plane and edge
%moment loading
Nx = input(' Write Nx value with N/mm =');
Ny = input(' Write Ny value with N/mm =');
Nxy = input(' Write Nxy value with N/mm =');
Mx = input(' Write Mx value with Nm/mm =');
My = input(' Write My value with Nm/mm =');
Mxy = input(' Write Mxy value with Nm/mm =');
N = [Nx;Ny;Nxy]
M = [Mx;My;Mxy]
% sigma_x= Nx / h;
% sigma_y= Ny / h;
% shear_xy= Nxy / h;
% stresses= [sigma_x;sigma_y;shear_xy] 
Q11= E_1 / (1-(v_12*v_21));
Q12= v_12 * E_2 / (1-(v_12*v_21));
Q22 = E_2 / (1-(v_12*v_21));
Q66 = G_12;
% zbar value for k=1,2,3,4 ...
for i = 1:Nplies;
  zbar(i) = (h+ h_ply)/2 - i*h_ply ;
end;

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
NT = zeros(3,1);
MT = zeros(3,1);
%Calculation Q , Qbar , A ,B and D matrices
for i = 1:Nplies;
  theta  = thetadb(i) * pi / 180; % lamina i angle in radians, from bottom
  m = cos(theta) ;
  n = sin(theta) ; 
Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66] ;
Q_11 = Q11 * (m^4) + 2 * ( Q12+ 2 * Q66)* m^2 * n^2 + Q22 * n^4 ;
Q_22 = Q11 * (n^4) + 2 * ( Q12+ 2 * Q66)* m^2 * n^2 + Q22 * m^4 ;
Q_12 = Q12 * (m^4) + Q12 * (n^4) + (Q11+ Q22 - 4*Q66) * m^2 * n^2;
Q_66 = (Q11+ Q22 - 2*Q12) * m^2 * n^2 + Q66 * (m^2 - n^2)^2 ;
Q_16 = -1 * m * n^3 * Q22 + m^3 *n * Q11 - m * n * (m^2 - n^2) * (Q12+ 2 * Q66);
Q_26 = -1 * m^3 * n * Q22 + m^3 *n * Q11 - m * n * (m^2 - n^2) * (Q12+ 2 * Q66);
format long
Qbar = [Q_11,Q_12,Q_16;Q_12,Q_22,Q_26;Q_16,Q_26,Q_66]; 
  Sbar = inv(Qbar);
  A = A + Qbar * h_ply;
  B = B + Qbar * h_ply * zbar(i); 
  D = D + Qbar * (h_ply * zbar(i)^2  + h_ply^3 / 12);
end;
ABD = [A,B;B,D] %stiffness matrix equation
abd = inv(ABD); %inverse of stiffness matrix
a = abd(1:3,1:3); 
b = abd(1:3,4:6);
d = abd(4:6,4:6);
%For strains
format long
load = [Nx;Ny;Nxy;Mx;My;Mxy] ; 
epsilonk =  abd * load 
epsilon0 = epsilonk(1:3)
K= epsilonk(4:6)

% for i = 1:Nplies;
% epsilont  = R +  zbar(i) .* H
% end
epsilont  = epsilon0 +  zbar(Nplies) * K %zbar(Nplies) for outer most fiber at top layer
epsilonb  = epsilon0 +  zbar(1) * K %zbar(Nplies) for outer most fiber at top layer
 T1 = [ m^2 n^2 m*n; n^2 m^2 -m*n; -m*n m*n (m^2 - n^2)]; %transformation matrix for strain
 T2 = [ m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n (m^2 - n^2)]; %transformation matrix for stress 
 %in plane stressess
epsilon1=epsilont';
strain= epsilon1 *T1
sigmas= Qbar .* epsilon1;
stresses= sigmas .* T2;
sigma1t= stresses(1,1)+stresses(1,2)+stresses(1,3);
sigma2t= stresses(2,1)+stresses(2,2)+stresses(2,3);
tau12t= stresses(3,1)+stresses(3,2)+stresses(3,3);
%for bottom layer
epsilon2=epsilonb';
strain2= epsilon2 * T1
sigmas2= Qbar .* epsilon2;
stresses2= sigmas2 .* T2;
sigma1b= stresses2(1,1)+stresses2(1,2)+stresses2(1,3);
sigma2b= stresses2(2,1)+stresses2(2,2)+stresses2(2,3);
tau12b= stresses2(3,1)+stresses2(3,2)+stresses2(3,3);
Qbar
Sbar
A % with unit N/mm
B % with unit N
D % with unit N*mm
a % with unit mm/N
b % with unit 1/N
d  % with unit 1/N*mm
sigma1t % with unit MPa
sigma2t % with unit MPa
tau12t  % with unit MPa
sigma1b % with unit MPa
sigma2b % with unit MPa
tau12b  % with unit MPa




%Arda Pamuk, Metu Ncc Aerospace Engineering Department 
