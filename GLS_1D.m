% a modified leap-frog scheme to solve Ginzburg-Landau-Schrodinger equation
function psin = GLS_1D(N,tau,gam)
% N = 64; 
le = 20;
h = le/N;
T = 1; 
x = -le/2:h:le/2; x = x(2:end-1)';
c = @(k,gam) ((-1)^k)*gamma(gam+1)/(gamma(gam/2-k+1)*gamma(gam/2+k+1));
alpha = 1; beta = 1;
deltax_gam = coef_matrix(N,gam,c,h);
init = sech(x).*exp(3i*x);
% calculate psi1
nonlinear_term = abs(init).^2 -1;
nonlinear_matrix = diag(nonlinear_term);
left_matrix = (alpha-beta*1i)/tau*eye(N-1) + deltax_gam/2 + nonlinear_matrix/2;
right_matrix = (alpha-beta*1i)/tau*eye(N-1) - deltax_gam/2 - nonlinear_matrix/2;
psi1 = left_matrix\(right_matrix*init);

psin_pre = init; psin = psi1; 
mass0 = h*sum(abs(init).^2);
energy0 = 2*h*real(sum((deltax_gam*init).*conj(init))) + h*real(sum(abs(abs(init).^2 -1).^2));
mass = [mass0];
energy = [energy0];
% calculate psin
num_steps = round(T/tau);
for n = 1:num_steps-1
    tn = n*tau;
    mass_psin = h*sum(abs(psin).^2);
    energy_psin = h*real(sum((deltax_gam*psin).*conj(psin))) + h*real(sum((deltax_gam*psin_pre).*conj(psin_pre)))+...
                    h*real(sum((abs(psin).^2 -1).*conj(abs(psin_pre).^2 -1)));
    mass = [mass mass_psin];
    energy = [energy energy_psin];
    nonlinear_term = abs(psin).^2 -1;
    nonlinear_matrix = diag(nonlinear_term);
    left_matrix = (alpha-beta*1i)/(2*tau)*eye(N-1) + deltax_gam/2 + nonlinear_matrix/2;
    right_matrix = (alpha-beta*1i)/(2*tau)*eye(N-1) - deltax_gam/2 - nonlinear_matrix/2;
    psin_next = left_matrix\(right_matrix*psin_pre);
    psin_pre = psin;
    psin = psin_next;
end
end











