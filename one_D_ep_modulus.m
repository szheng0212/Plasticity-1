clear;
n_max=100; %numbers of ciculation
%initialization of the strain variables
epsilon=zeros(n_max,1);
epsilon_p=zeros(n_max,1);
delta_eps=0.0025; 
%initialization of the stress variables
sigma=zeros(n_max,1);
sigma_trial=zeros(n_max,1);
sigma_y=1.0e5; 
alpha=zeros(n_max,1);
alpha_trial=zeros(n_max,1);
f_trial=zeros(n_max+1,1);
E=1.0e6;K=0.25e6; % Young's modulus
for n=1:n_max
    epsilon(n+1)=epsilon(n)+delta_eps;
    sigma_trial(n+1)=E*(epsilon(n+1)-epsilon_p(n));
    alpha_trial(n+1)=alpha(n);
    f_trial(n+1)=abs(sigma_trial(n+1))-sigma_y-K*alpha_trial(n+1);
    if f_trial(n+1)<=0
        sigma(n+1)=sigma_trial(n+1);
        alpha(n+1)=alpha_trial(n+1);
    else
        sigma(n+1)=sigma(n)+(E*K)/(E+K)*delta_eps;
        alpha(n+1)=alpha(n)+f_trial(n+1)/(E+K);
    end        
end
figure (1)
plot(epsilon,sigma)
figure(2)
plot(epsilon,alpha)