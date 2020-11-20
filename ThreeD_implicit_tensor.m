clc;clear;
%Define constants
global G K H c0 eta eta_bar xi delta_gamma;
G=7.69e6; v=0.3;
K=2*G*(1+v)/(3-6*v);
H=-100e3;
n_max=100;
inc_e_1=-0.0001;
stress_c=-50e3;          %Confining pressure unit:Pa
c0=10e3;
e_tol=1.0e-3;%tolerant error
phi=0/180*pi;%frictional angle
eta=6*sin(phi)/sqrt(3)/(3+sin(phi));
xi=6*cos(phi)/sqrt(3)/(3+sin(phi));
psi=0/180*pi;%dilantancy angle
eta_bar=6*sin(psi)/sqrt(3)/(3+sin(psi));

path=3;    % 1 for oedometer test; 2 for undrained triaxial compression; 3 for simple shear
delta_eps=zeros(3,3);
if path==1         % oedometer test
    delta_eps(1,1)=inc_e_1*1.0;
elseif path==2     % undrained triaxial compression
    delta_eps(1,1)=inc_e_1*1.0;
    delta_eps(2,2)=-inc_e_1*0.5;
    delta_eps(3,3)=-inc_e_1*0.5;
elseif path==3     % simple shear
    delta_eps(1,3)=inc_e_1*0.5;
    delta_eps(3,1)=inc_e_1*0.5;
end

% Define the strain arrays
epsilon=zeros(3,3,n_max);
epsilon_e=zeros(3,3,n_max);
epsilon_p_bar=zeros(n_max,1);
epsilon_e0=zeros(3,3);

%Define the stress arrays
stress=zeros(3,3,n_max);
stress0=eye(3)*stress_c;
stress(:,:,1)=stress0;
p=zeros(n_max,1);
p(1)=trace(stress0)/3;
s=zeros(3,3,n_max);
s(:,:,1)=get_dev_tensor(stress0);
sqrt_J2=zeros(n_max,1);

for n=1:n_max
    delta_gamma=0;% initialization of delta_gamma
    epsilon_e_trial=epsilon_e(:,:,n)+delta_eps;
    epsilon(:,:,n+1)=epsilon(:,:,n)+delta_eps;
    epsilon_p_bar_trial=epsilon_p_bar(n);
    s_trial=s(:,:,1)+2*G*get_dev_tensor(epsilon_e_trial);
    p_trial=p(1)+K*trace(epsilon_e_trial);
    phi_wave=get_residual(s_trial,p_trial,epsilon_p_bar_trial);
    if phi_wave<=0
        s(:,:,n+1)=s_trial;
        p(n+1)=p_trial;
        stress(:,:,n+1)=s_trial+p_trial*eye(3);
        sqrt_J2(n+1)=sqrt(0.5*double_dot(s_trial,s_trial));
        epsilon_e(:,:,n+1)=epsilon_e_trial;
        epsilon_p_bar(n+1)=epsilon_p_bar_trial;       
    else
        while 1
            d=-G-K*eta*eta_bar-xi*xi*H;
            delta_gamma=delta_gamma-phi_wave/d;
            epsilon_p_bar(n+1)=epsilon_p_bar(n)+xi*delta_gamma;
            phi_wave=get_residual(s_trial,p_trial,epsilon_p_bar(n+1));
            if abs(phi_wave)<=e_tol
                s(:,:,n+1)=(1-G*delta_gamma/sqrt(0.5*double_dot(s_trial,s_trial)))*s_trial;
                p(n+1)=p_trial-K*eta_bar*delta_gamma;
                stress(:,:,n+1)=s(:,:,n+1)+ p(n+1)*eye(3);
                epsilon_e(:,:,n+1)=epsilon_e0+(s(:,:,n+1)-s(:,:,1))/(2*G)+(p(n+1)-p(1))/(3*K)*eye(3);
                sqrt_J2(n+1)=sqrt(0.5*double_dot(s(:,:,n+1),s(:,:,n+1)));
                break
            end                                 
        end
    end        
end

%post-process
t=1:n_max;
if path==1
    plot_oedometer_test_tensor
elseif path==2
    plot_undrained_triaxial_compression_tensor
elseif path==3
    plot_simple_shear_tensor
end


function dev_tensor=get_dev_tensor(tensor)
    dev_tensor=tensor-eye(3)*trace(tensor)/3;
end

function phi_wave=get_residual(s,p,e_p_bar)
    global delta_gamma G K eta eta_bar xi;
    J2=0.5*double_dot(s,s);
    c=get_c(e_p_bar);
    phi_wave=sqrt(J2)-G*delta_gamma+eta*(p-K*eta_bar*delta_gamma)-xi*c;
end

function c=get_c(e_p_bar)
    global H c0;
    c=c0+H*e_p_bar;
end

function output=double_dot(tensor1,tensor2)
    output=sum(dot(tensor1,tensor2));
end
