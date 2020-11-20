clc;clear;
%Define constants
global n_max G K H c0 eta eta_bar xi delta_gamma;
G=7.69e6; v=0.3;
K=2*G*(1+v)/(3-6*v);
H=0e3;
n_max=100;
inc_e_1=-0.0001;
stress_c=-50e3;       %Confining pressure unit:Pa
c0=10e3;
e_tol=1.0e-3;%tolerant error
phi=30/180*pi;%frictional angle
eta=6*sin(phi)/sqrt(3)/(3+sin(phi));
xi=6*cos(phi)/sqrt(3)/(3+sin(phi));
psi=0/180*pi;%dilantancy angle
eta_bar=6*sin(psi)/sqrt(3)/(3+sin(psi));

path=1;    % 1 for oedometer test; 2 for undrained triaxial compression; 3 for simple shear
% Define the strain arrays
if path==1
    delta_eps=[inc_e_1*1.0,0.0,0.0,0.0,0.0,0];% oedometer test
elseif path==2
    delta_eps=[inc_e_1*1.0,-inc_e_1*0.5,-inc_e_1*0.5,0.0,0.0,0];% undrained triaxial compression
elseif path==3
    delta_eps=[0.0,0.0,0.0,inc_e_1*1.0,0.0,0];% simple shear
end

epsilon=zeros(n_max,6);
epsilon_e=zeros(n_max,6);
epsilon_p_bar=zeros(n_max,1);
epsilon_e0=[0,0,0,0,0,0];
epsilon_e(1,:)=epsilon_e0;

%Define the stress arrays
stress=zeros(n_max,6);
stress0=[stress_c,stress_c,stress_c,0,0,0];
stress(1,:)=stress0;
s=zeros(n_max,6);
s0=get_dev_stress(stress0);
s(1,:)=s0;
p=zeros(n_max,1);
p0=get_mean_stress(stress0);
p(1)=p0;
sqrt_J2=zeros(n_max,1);

for n=1:n_max
    delta_gamma=0;% initialization of delta_gamma
    epsilon_e_trial=epsilon_e(n,:)+delta_eps;
    epsilon(n+1,:)=epsilon(n,:)+delta_eps;
    epsilon_p_bar_trial=epsilon_p_bar(n);
    s_trial=s0+2*G*get_dev_strain(epsilon_e_trial);
    p_trial=p0+K*get_vol_strain(epsilon_e_trial);
    phi_wave=get_residual(s_trial,p_trial,epsilon_p_bar_trial);
    if phi_wave<=0
        s(n+1,:)=s_trial;
        p(n+1)=p_trial;
        stress(n+1,:)=s_trial+[p_trial,p_trial,p_trial,0,0,0];
        sqrt_J2(n+1)=sqrt(get_J2(s(n+1,:)));
        epsilon_e(n+1,:)=epsilon_e_trial;
        epsilon_p_bar(n+1)=epsilon_p_bar_trial;       
    else
        delta_gamma=(eta*K*trace(delta_eps)+(G/sqrt(get_J2(s(n,:))))*double_dot(s(n,:),delta_eps))/(K*eta*eta_bar+G);
        delta_s=2*G*get_dev_strain(delta_eps)-delta_gamma*G/sqrt(get_J2(s(n,:)))*s(n,:);
        s(n+1,:)= s(n,:)+delta_s;
        delta_p=K*trace(delta_eps)-delta_gamma*K*eta_bar;
        p(n+1)=p(n)+delta_p;
        stress(n+1,:)=s(n+1,:)+[p(n+1),p(n+1),p(n+1),0,0,0];
        epsilon_e(n+1,:)=epsilon_e0+(s(n+1,:)-s0)/(2*G)+(p(n+1)-p0)/(3*K)*[1,1,1,0,0,0];
        sqrt_J2(n+1)=sqrt(get_J2(s(n+1,:)));
    end        
end

%post-process
t=1:n_max;
if path==1
    plot_oedometer
elseif path==2
    plot_undrained_triaxial_compression
elseif path==3
    plot_simple_shear
end

function dev_stress=get_dev_stress(stress)
p=get_mean_stress(stress);
dev_stress=stress;
dev_stress(1)=stress(1)-p;
dev_stress(2)=stress(2)-p;
dev_stress(3)=stress(3)-p;
end

function dev_strain=get_dev_strain(strain)
    dev_strain=zeros(1,6);
    trace=strain(1)+strain(2)+strain(3);
    dev_strain(1)=strain(1)-trace/3;
    dev_strain(2)=strain(2)-trace/3;
    dev_strain(3)=strain(3)-trace/3;
    dev_strain(4)=strain(4);
    dev_strain(5)=strain(5);
    dev_strain(6)=strain(6);
end

function p=get_mean_stress(stress)
p=(stress(1)+stress(2)+stress(3))/3;
end

function vol_strain=get_vol_strain(strain)
    vol_strain=strain(1)+strain(2)+strain(3);
end

function phi_wave=get_residual(s,p,e_p_bar)
    global delta_gamma G K eta eta_bar xi;
    J2=get_J2(s);
    c=get_c(e_p_bar);
    phi_wave=sqrt(J2)-G*delta_gamma+eta*(p-K*eta_bar*delta_gamma)-xi*c;
end

function J2=get_J2(s)
    J2=0.5*(s(1)*s(1)+s(2)*s(2)+s(3)*s(3)+2*(s(4)*s(4)+s(5)*s(5)+s(6)*s(6)));
end

function c=get_c(e_p_bar)
    global H c0;
    c=c0+H*e_p_bar;
end

function tr=trace(tensor)
    tr=tensor(1)+tensor(2)+tensor(3);
end
function prod=double_dot(tensor1,tensor2)
    prod=tensor1(1)*tensor2(1)+tensor1(2)*tensor2(2)+tensor1(3)*tensor2(3)+2*(tensor1(4)*tensor2(4)+tensor1(5)*tensor2(5)+tensor1(6)*tensor2(6));
end