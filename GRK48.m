% calculate reference solution for 1D Ginzburg-Landau-Schrodinger equation
% use the fixed-point iteration


tau = 1e-2;
le = 20; 
N = 256;
h = le/N;
T = 1; 
x = -le/2:h:le/2; x = x(2:end-1)';
c = @(k,gam) ((-1)^k)*gamma(gam+1)/(gamma(gam/2-k+1)*gamma(gam/2+k+1));
gam = 1.8;
alpha = 1; beta = 1;
deltax_gam = coef_matrix(N,gam,c,h);
init = sech(x).*exp(3i*x);

w1 = 1/8 - sqrt(30)/144; 
w2 = 0.5*(sqrt((15+2*(sqrt(30)))/35));
w3 = w2*(1/6 + sqrt(30)/24);
w4 = w2*(1/21 + 5*sqrt(30)/168);
w5 = w2 - 2*w3;
w11 = 1/8 + sqrt(30)/144;
w22 = 0.5*(sqrt((15-2*(sqrt(30)))/35));
w33 = w22*(1/6 - sqrt(30)/24);
w44 = w22*(1/21 - 5*sqrt(30)/168);
w55 = w22 - 2*w33;

A = [w1 w11-w3+w44 w11-w3-w44 w1-w5
     w1-w33+w4 w11 w11-w55 w1-w33-w4
     w1+w33+w4 w11+w55 w11 w1+w33-w4
     w1+w5 w11+w3+w44 w11+w3-w44 w1];
B = [2*w1 2*w11 2*w11 2*w1];

num_steps = round(T/tau);
Lf = -deltax_gam/(alpha - beta*1i); A_kron_Lf = kron(A,Lf);
matrix_left = speye(4*(N-1)) - tau*A_kron_Lf;
Un = init;
for n = 1:num_steps
    iter_err = 1;  iter_num = 0;
    Uni_pre = [Un,Un,Un,Un]; % 四列
    Un_store = [Un;Un;Un;Un]; % 竖着放四�?
    while ((iter_err > 1e-14) && (iter_num < 20))
        Fni_store = ((1-abs(Uni_pre).^2).*Uni_pre)/(alpha - beta*1i);
        AF = Fni_store*A';
        AF = AF(:);
        Uni_pre = Uni_pre(:);
        Uni_next = matrix_left\(Un_store + tau*AF);
        iter_err = max(abs(Uni_next-Uni_pre));
        iter_num = iter_num + 1;
        Uni_pre = Uni_next;
        Uni_pre = reshape(Uni_pre,N-1,4);
    end
    Uni_store = Uni_next;
    Uni_store = reshape(Uni_store,N-1,4);
    Fni_store = ((1-abs(Uni_store).^2).*Uni_store)/(alpha - beta*1i);
    Un_next = Un + tau*Lf*Uni_store*(B') + tau*Fni_store*(B');
    Un = Un_next;    
    tn = n*tau
end
UR = Un;
% save 1c2.mat UR
% save 2.mat UR
plot(x,abs(UR))

