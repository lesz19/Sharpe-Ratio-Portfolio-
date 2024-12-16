clear all;
clc;
close all;

options = optimset('Display','iter', 'FunValCheck', 'off', 'Algorithm', 'active-set');

%Estados
N = 6;   
%Acciones
M = 3;  


pi1(:,:,1) = [    0.8147    0.2785    0.9572    0.7922    0.6787    0.0060;
    0.9058    0.7469    0.9854    0.5595    0.7577    0.0318;
    0.1270    0.9575    0.0003    0.6557    0.7431    0.2769;
    0.0134    0.8649    0.1419    0.0357    0.3922    0.0462;
    0.6324    0.1576    0.4218    0.0491    0.0555    0.0971;
    0.9975    0.9706    0.9157    0.7340    0.1712    0.0235];   
    
pi1(:,:,2) = [    0.0948    0.0655    0.0094    0.1190    0.0513 0.0472;
    0.6171    0.9952    0.9547    0.0984    0.2551    0.0386;
    0.4502    0.9869    0.3760    0.0597    0.0060    0.1493;
    0.7100    0.5998    0.0797    0.0404    0.0091    0.2575; 
    0.7387    0.8456    0.3551    0.1853    0.0909    0.0407;
    0.8816    0.6883    0.2106    0.0238    0.0593    0.1543]; 

pi1(:,:,3) = [    0.8143    0.6160    0.9172    0.0759    0.5688    0.3112;
    0.6435    0.8733    0.2858    0.0540    0.4694    0.0285;
    0.9293    0.3517    0.7572    0.5308    0.0119    0.1656;
    0.3500    0.0308    0.0537    0.7792    0.3371    0.6020;
    0.1966    0.5853    0.3804    0.0340    0.1622    0.2630;
    0.2511    0.4497    0.5678    0.1299    0.0943    0.0541];


u1(:,:,1) = [41    14    48    40    34    36;
    46    28    25    48    38     2;
     7    48    41    33    38    14;
    46    49     8     2    12     3;
    32     8    22    43    33     5;
     5    49    46    47     9    42];

u1(:,:,2) = [    90    29    18    21    75    14;
    12    42    36    91    74    22;
    99    47     6    68    57    90;
    54    77    53    47    19     8;
    71    82    34    92    60    25;
   100    11    18    11    30     6];

u1(:,:,3) = [45    46    30     9    34    61;
     2    11     5    78    30    53;
    90   100    51    91    75    73;
    20    34    77    54     2    71;
    10    30    64    11     5    79;
    31     7     9    83    67    29];

%Cálculo de las matrices W para Método
W1=zeros(N,M);

for i1=1:N
        for k1=1:M
                for j1=1:M
                         W1(i1,k1)=W1(i1,k1)+u1(i1,j1,k1)*pi1(i1,j1,k1);
                end
        end
end

NM =  N*M;

WW1=zeros(NM,1);
j1 = 0;
for i1=1:N
         for k1=1:M
                 j1 = j1 + 1;
                 WW1(j1)=W1(i1,k1);
         end
end


    A = zeros(N, NM);

    for j = 1 : N
        aux = 0;
        for k = 1 : M
            for i = 1 : N
                aux = aux + 1;            
                if i == j
                    A(j,aux) = pi1(i,j,k) - 1;
                else
                    A(j,aux) = pi1(i,j,k);
                end
            end
        end
    end    

    a1col = size(A,2);
    A_aux = ones(1,a1col);

    Aeq = [A;A_aux];
    beq = [zeros(N,1);1]; 

    x1 = zeros(1,NM);    
    x = ones(NM,1)/(NM);     

nf = 10000;%550000; %5000
n0 = 3100;%100

gamma0 = 3e-9;
gamma0u = 0.09;
gamma0v = 9e-7;

% teta0 = 1e-9;
delta0 = 1e-2;
v00 = 0.7;
v10 = ones(NM,1)*.7;
mu00 = ones(NM,1)*.7;
mu10 = ones(N,1)*.7; % listo

gamma = gamma0;
gammau = gamma0u;
gammav = gamma0v;

delta = delta0;
v0 = v00;
v1  = v10;
mu0 = mu00;
mu1 = mu10;
alpha = 10;

X1 = zeros(NM,nf);
Ng = zeros(1,nf);
P = zeros(1,nf);
D = zeros(1,nf);
R = zeros(1,nf);
F = zeros(1,nf);

X1(:,1) = x;


    f1 = 0;
    f2 = 0;
    for j1=1:NM
        f1 = f1 + WW1(j1)*x(j1);
        f2 = f2 + (WW1(j1)^2*x(j1)-WW1(j1)^2*x(j1)^2);
    end   
    
    lambda = 0.01;        
%     P(1) = 0;
%     D(1) = f2;
%     R(1) = f1;
%     F(1) = f1 - lambda*f2 + delta*sum(x)^2;

    
for n = 1 : nf   

     if n > n0
         gamma = gamma0/(1+n-n0)^(1/2);
         gammau = gamma0u/(1+n-n0)^(1/2);  
         gammav = gamma0v/(1+n-n0)^(1/2);
         delta = delta0*(1+log(n-n0))/(1+n-n0)^(1/2);
     end 

 % Gradient of theta by x 
    df1 = zeros(NM,1); 
    df2 = zeros(NM,1); 
        for j = 1 : NM
            df1(j) =  df1(j) +WW1(j1);
            df2(j) =  df2(j) + (WW1(j1)^2*x(j1)-WW1(j1)^2*x(j1)^2)^(-1/2)*(WW1(j)^2-WW1(j)^2*x(j));              
        end 
    e = ones(NM,1);
    G = v0 * df2 - v0 * alpha * df1 - v1 + mu0 + A'*mu1-delta*x;    
    x = x - gamma*G;  
    
 % Gradient of alpha     
    dfa = 0;
        for j = 1 : NM
            dfa =  dfa + WW1(j)*x(j);
        end
    Ga = 1 + v0 * dfa + delta*alpha;
    alpha = alpha - gamma*Ga;  
    alpha = max(0,alpha);

% Gradient of v0  
    df1 = zeros(NM,1); 
    df2 = zeros(NM,1); 
        for j = 1 : NM
            df1(j) =  df1(j) +WW1(j1)*x(j1);
            df2(j) =  df2(j) + (WW1(j1)^2*x(j1)-WW1(j1)^2*x(j1)^2)^(1/2);              
        end 
    df = df2-alpha*df1;
    dfv0 = sum(df);       
    Gv0 = dfv0 - delta*v0;
    v0 = v0 - gammav*Gv0;
    v0 = max(0,v0);

% Gradient of v1    
    dfv1 = sum(x); 
    Gv1 = -dfv1*e + delta*v1;
    v1 = v1 - gammav*Gv1;
    v1 = max(0,v1);

% Gradient of mu0     
    dfmu0 = zeros(NM,1);
    Gmu0 = (sum(x)-1)*e - delta*mu0;
    mu0 = mu0 - gammau*Gmu0;
    mu0 = max(0,mu0);

% Gradient of mu1     
    dfmu1 = zeros(N,1);
    Gmu1 =  A*x - delta*mu1;
    mu1 = mu1 - gammau*Gmu1;
    mu1 = max(0,mu1);
       

    Ng(n) = norm(Aeq*x-beq);
    f1 = 0;
    f2 = 0;
    for j1=1:NM
        f1 = f1 + WW1(j1)*x(j1);
        f2 = f2 + (WW1(j1)^2*x(j1)-WW1(j1)^2*x(j1)^2);
    end    
    
    X1(:,n) = x;  
    P(n) = f2/f1;   

    R(n) = f1;
    D(n) = f2^(1/2);
    F(n) = f1 - lambda*f2 + delta*sum(x)^2;

end

figure 
hold on
plot(F(:),'-k','LineWidth',2)
xlabel('Iterations')
title('Functional')
grid on
legend('E[R(\omega)]-Var(R(\omega))')
hold off

 t=1:nf; 
 
figure;
hold on
plot(t,X1(1,:),'--k','LineWidth',2)
plot(t,X1(2,:),'--b','LineWidth',2)
plot(t,X1(3,:),':y','LineWidth',2)
plot(t,X1(4,:),'-.m','LineWidth',2)
plot(t,X1(5,:),':g','LineWidth',2)
plot(t,X1(6,:),'-.c','LineWidth',2)
plot(t,X1(7,:),'--r','LineWidth',2)
plot(t,X1(8,:),'--b','LineWidth',2)
plot(t,X1(9,:),'-.r','LineWidth',2)
plot(t,X1(10,:),'-k','LineWidth',2)
plot(t,X1(11,:),'-.b','LineWidth',2)
plot(t,X1(12,:),'-y','LineWidth',2)
plot(t,X1(13,:),':m','LineWidth',2)
plot(t,X1(14,:),'--g','LineWidth',2)
plot(t,X1(15,:),'--c','LineWidth',2)
plot(t,X1(16,:),'-.r','LineWidth',2)
plot(t,X1(17,:),'-.k','LineWidth',2)
plot(t,X1(18,:),'--r','LineWidth',2)
xlabel('Iterations')
ylabel('Probability')
title('Strategies \theta')
hold on
legend('\theta(1,1)','\theta(2,1)','\theta(3,1)','\theta(4,1)','\theta(5,1)','\theta(6,1)','\theta(1,2)','\theta(2,2)','\theta(3,2)','\theta(4,2)','\theta(5,2)','\theta(6,2)','\theta(1,3)','\theta(2,3)','\theta(3,3)','\theta(4,3)','\theta(5,3)','\theta(6,3)');
hold off



figure;
count2=timeseries(Ng(:),t);
plot(count2,'-r','MarkerSize',16,'LineWidth',2)
title('Norm of Constraints')
hold on
ylabel('')
% hold on
% h = legend('Constraints')

figure;
count2=timeseries(P(:),t);
plot(count2,'-b','MarkerSize',16,'LineWidth',2)
title('Sharpe-ratio portfolio')
grid on
xlabel('Iterations')
ylabel('')
legend('Var(R(\omega))^{1/2}/E[R(\omega)]')
hold off

figure
hold on
plot(D(:),R(:),'-k','LineWidth',2)
xlabel('Standard deviation')
ylabel('Utility')
title('Eficient Frontier')
grid on
legend('Var(R(\omega))^{1/2} vs E[R(\omega)]')
hold off

norm(Aeq*x-beq)
cik = zeros(N,M);
j = 0;
for i = 1 : N
    for k = 1 : M
        j = j + 1;
        cik(i,k) = x(j);        
    end

end

p = zeros(1,N);
for i = 1 : N
    for k = 1 : M
        p(i) = p(i) + cik(i,k);         
    end
end
p
sum(p)

% varaible \pi(a|s) por dik
dik = zeros(N,M);
for i = 1 : N
    for k = 1 : M
        dik(i,k) = cik(i,k)/p(i);         
    end
end
dik
