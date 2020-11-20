clear;
n_max=10; %numbers of ciculation
%initialization of the strain variables
epsilon=zeros(n_max,1);
epsilon_p=zeros(n_max,1);
delta_eps=0.025; 
%initialization of the stress variables
sigma=zeros(n_max,1);
sigma_trial=zeros(n_max,1);
sigma_y=1.0e5; 
alpha=zeros(n_max,1);
alpha_trial=zeros(n_max,1);
f_trial=zeros(n_max+1,1);
E=1.0e6;K=0e6; % Young's modulus
for n=1:n_max
    epsilon(n+1)=epsilon(n)+delta_eps;
    sigma_trial(n+1)=E*(epsilon(n+1)-epsilon_p(n));
    alpha_trial(n+1)=alpha(n);
    f_trial(n+1)=abs(sigma_trial(n+1))-sigma_y-K*alpha_trial(n+1);
    if f_trial(n+1)<=0
        sigma(n+1)=sigma_trial(n+1);
    else
        delta_gamma=f_trial(n+1)/(E+K);
        sigma(n+1)=(1-delta_gamma*E/abs(sigma_trial(n+1)))*sigma_trial(n+1);
        epsilon_p(n+1)=epsilon_p(n)+delta_gamma*sigma_trial(n+1)/abs(sigma_trial(n+1));
        alpha(n+1)=alpha(n)+delta_gamma;
    end        
end
plot(epsilon,sigma)