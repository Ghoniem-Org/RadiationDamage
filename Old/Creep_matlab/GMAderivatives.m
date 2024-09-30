function F = GMAderivatives(t,x)

% x is a 4 slot vector [rho_m, rho_s, rho_b, R_sb]
rho_m = x(1);                                  %[#of dislocations/m^2]density of mobile dislocations
rho_s = x(2);                                  %[#of dislocations/m^2]density of static dislocations
rho_b = x(3);                                  %[#of dislocations/m^2]density of boundary dislocations
R_sb = x(4);                                   %[m]subgrain radius
tau_oct = 150e6;                               %[Pa]applied stress
%---------- Parameters -------------%
b = 2e-10;%1e-10;                                      %[m]burgers vector
D_o = 2e-5;                                     %[m^2/s] Diffusion coefficient
SFE = 0.2;                                      %[J/m^2]Stacking fault energy
delta = 4e-9;                                   %[m]seperation between dislocations for spontaneous annihilation
sigma_o = 20e6;                                 %[Pa]friction stress    
E = 1.77e11;                                    %[Pa]Young's modulus
nu = 0.3;                                       %Poisson's ratio 
N_p =1.63e16; % 1.63e22;                         %[#precipitates/m^3] volume concentration of precipitates
r_p = 405e-10;                                  %[m] precipitates mean radius 
k = 1.38065e-23;                                %[m^2*kg/(s^2*K)]or[J/K]Boltzmann constant 
Omega = 1.19e-29;                               %[m^3] atomic volume
T = 550 + 273;                                  %[K] Temperature
mu = E/(2*(1 + nu));                            %[Pa] Shear modulus 
W_g = 6e-19;                                    %Fitting parameter for glide velocity (eq12 on 1990-ghon-matt paper)
a1 = 1e11;                                      %Fitting parameter for glide velocity (eq12 on 1990-ghon-matt paper)
c_jog = 0.0389;                                 %[?]jog concentration
eta_v = 10^3*1*c_jog*b*(SFE/(mu*b))^2;          %[?] parameter defining the transfer of defects to jog on dislocations (eq23 on 1990-ghon-matt paper)
E_core = 1.3*1.602e-19;                         %[J]
E_s = 2.8*1.602e-19;                            %[J]
E_m = E_s/2;
K_c = 10;                         %constant from Holt analysis (eq44 on 1990-ghon-matt paper)
D_p = D_o*exp(-E_core/(k*T));     %[m^2/s]core diffusion coefficient
D_v = D_o*exp(-E_m/(k*T));        %[m^2/s]vacancy diffusion coefficient
D_s = D_o*exp(-E_s/(k*T));        %[m^2/s]Lattice self diffusion coefficient 
Beta = 1.15e5;                    %Fitting varialbe controlling density of sources 
xi = 1;                           %Uniltless constant (eq7 on 1990-ghon-matt paper)
zeta = 0.425;                     %constant (eq49 on 1990-ghon-matt paper)
%-----------------------------------%
h = 1/(R_sb*(rho_s + rho_b));                          %[m]dislocations spacing within the subgrain walls(eq3 on 1990-ghon-matt paper)
gamma_sb = mu*b^2*rho_b*R_sb/3;                        %[J/m^2]Low-angle subgrain boundry energy per unit area(eq36 on 1990-ghon-matt paper)
p_s = (4/3)*mu*b^2*rho_b;                              %[Pa]Pressure for subgrain growth(Gibbs-Thompson effect)(eq38 on 1990-ghon-matt paper)
M_sb_c = 2*pi*b*D_p*Omega/(h^2*k*T);                   %unit?[m^3/(N.s)]Core mobility(eq39 on 1990-ghon-matt paper)
M_sb_L = 2*pi*eta_v*D_v*Omega/(b*k*T);                 %unit?[m^3/(N.s)]Lattice mobility(eq40 on 1990-ghon-matt paper)         
M_sb = M_sb_c + M_sb_L;                                %unir?[m^3/(N.s)]
lambda_d = 1/sqrt(rho_m);                              %[m]inter-dislocation spacing
lambda_p = 1/sqrt(N_p*r_p);                            %[m]inter-precipitates spacing
lambda = 1/(lambda_d^-1 + lambda_p^-1);                %[m]inter-obstacles spacing considering both dislocations and precipitates
alpha = (1/(pi*(1-nu)));                               
sigma_i = mu*b/(2*pi*lambda) + xi*mu*b*sqrt(rho_s);                        %[Pa]long range internal stress(eq7 on 1990-ghon-matt paper) 
sigma_e = tau_oct - sigma_i - sigma_o;                                     %[Pa]Effective stress on dislocations
sigma_s = mu*b*sqrt(rho_s)*alpha;                                          %[Pa]Static Dislocations stress???
sigma_m = mu*b*sqrt(rho_m)*alpha;                                          %[Pa]Mobile Dislocations stress???
v_g = (sigma_e*Omega*a1/(k*T))*exp(-W_g/(k*T));                            %[?]Empirical glide velocity (eq12 on 1990-ghon-matt paper)                      
v_cs = c_jog*2*pi*(D_s/b)*(Omega*sigma_s/(k*T))/log(1/sqrt(b/R_sb));       %[?]
v_cm = c_jog*2*pi*(D_s/b)*(Omega*sigma_m/(k*T))/log(1/sqrt(b^2*rho_m));    %[?]
%-----------------------------------%
R1 = sqrt(rho_m)^3;
R2 = Beta/(h^2*R_sb);
R3 = rho_m/(2*R_sb);
R4 = 8*sqrt(rho_m)^3*v_cm/v_g;
R5 = delta*rho_m*(rho_m + rho_s);
R6 = 8*(rho_s/h)*(v_cs/v_g);
R7 = delta*rho_m*rho_s;
R8 = 8*(1 - 2*zeta)*rho_s*(v_cs/h);
R9 = (rho_b/R_sb)*M_sb*(p_s - 2*pi*r_p^2*N_p*gamma_sb);
R10 = M_sb*(p_s - 2*pi*r_p^2*N_p*gamma_sb);
R11 = mu*eta_v*K_c*R_sb*(sqrt(rho_m + rho_s) - K_c/(2*R_sb))*(Omega*D_s/(k*T));
%---------- GMA eqns ---------------%
f1 = v_g*(R1 + R2 - R3 - R4 - R5);       % rho_m_dot
f2 = v_g*(R3 - R6 - R7);                 % rho_s_dot
f3 = R8 - R9;                            % rho_b_dot
f4 = R10 - R11;                          % R_sb_dot
f5=rho_m*b*v_g;
F = [f1; f2; f3; f4; f5]; 


    
    
    