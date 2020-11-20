% Define the strain arrays
n_max=50;
delta_eps=[-0.01,0.01,0,0,0];
epsilon=zeros(n_max,5);
epsilon_e=zeros(n_max,5);
epsilon_p_bar=zeros(n_max,1);

%Define the stress arrays
stress=zeros(n_max,5);
s=zeros(n_max,5);
p=zeros(n_max,1);
sqrt_J2=zeros(n_max,1);

%Define constants
global G K H c0 eta eta_bar xi delta_gamma;
G=1e6; v=0.25;
K=2*G*(1+v)/(3-6*v);
H=0e5;
c0=1e5;
e_tol=1.0e-3;

phi=0/180*pi;%frictional angle
eta=3*tan(phi)/sqrt(9+12*tan(phi)^2);
xi=3/sqrt(9+12*tan(phi)^2);

psi=0/180*pi;%dilantancy angle
eta_bar=3*tan(psi)/sqrt(9+12*tan(psi)^2);

for n=1:n_max
    delta_gamma=0;% initialization of delta_gamma
    epsilon_e_trial=epsilon_e(n,:)+delta_eps;
    epsilon(n+1,:)=epsilon(n,:)+delta_eps;
    epsilon_p_bar_trial=epsilon_p_bar(n);
    s_trial=2*G*get_dev_strain(epsilon_e_trial);
    p_trial=K*get_vol_strain(epsilon_e_trial);
    phi_wave=get_residual(s_trial,p_trial,epsilon_p_bar_trial);
    if phi_wave<=0
        stress(n+1,:)=s_trial+[p_trial,p_trial,p_trial,0,0];
        epsilon_e(n+1,:)=epsilon_e_trial;
        epsilon_p_bar(n+1)=epsilon_p_bar_trial;
        s(n+1,:)=s_trial;
        p(n+1)=p_trial;
        sqrt_J2(n+1)=sqrt(get_J2(s(n+1,:)));
    else
        while 1
            d=-G-K*eta*eta_bar-xi*xi*H;
            delta_gamma=delta_gamma-phi_wave/d;
            epsilon_p_bar(n+1)=epsilon_p_bar(n)+xi*delta_gamma;
            phi_wave=get_residual(s_trial,p_trial,epsilon_p_bar(n+1));
            if abs(phi_wave)<=e_tol
                s(n+1,:)=(1-G*delta_gamma/sqrt(get_J2(s_trial)))*s_trial;
                p(n+1)=p_trial-K*eta_bar*delta_gamma;
                stress(n+1,:)=s(n+1,:)+[ p(n+1), p(n+1), p(n+1),0,0];
                epsilon_e(n+1,:)=s(n+1,:)/(2*G)+p(n+1)/(3*K)*[1,1,1,0,0];
                sqrt_J2(n+1)=sqrt(get_J2(s(n+1,:)));
                break
            end                                 
        end
    end        
end
t=1:n_max;
figure(1)
plot(t,stress(t,1),t,sqrt_J2(t))
figure(2)
plot(t,stress(t,1),t,stress(t,2))

%plot(epsilon(t,4),sqrt_J2(t))
function dev_strain=get_dev_strain(strain)
    dev_strain=zeros(1,5);
    trace=strain(1)+strain(2)+strain(3);
    dev_strain(1)=strain(1)-trace/3;
    dev_strain(2)=strain(2)-trace/3;
    dev_strain(3)=strain(3)-trace/3;
    dev_strain(4)=strain(4);
    dev_strain(5)=strain(5);
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
    J2=0.5*(s(1)*s(1)+s(2)*s(2)+s(3)*s(3)+s(4)*s(4)+s(5)*s(5));
end

function c=get_c(e_p_bar)
    global H c0;
    c=c0+H*e_p_bar;
end