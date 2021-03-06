global P m c0 H eta xi eta_bar K De;
a=[0.25,0.25,0.5]; % material parameters
b=[1/3,1/3,2/3];
P1=[2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3]; P2=2*[a(1) a(2); a(2) a(1)]; P3=2*a(3)*eye(2);
Q1=[2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3]; Q2=1.5*[b(1) b(2); b(2) b(1)];Q3=1.5*b(3)*eye(2);
P=blkdiag(P1,P2,P3);
Q=blkdiag(Q1,Q2,Q3);
m=zeros(7,1); m(1)=1/3; m(2)=1/3; m(3)=1/3;

G=10e6;Gc=0.5*G;nu=0.25;K=2*G*(1+nu)/(3-6*nu); lamda=2*G*nu/(1-2*nu); %shear and cosserat modulus as well as poisson's ratio and lamda
De1=lamda*ones(3,3)+2*G*eye(3);  De2=[G+Gc G-Gc;G-Gc G+Gc];  De3=2*G*eye(2);
De=blkdiag(De1,De2,De3);
c0=10e3;
H=-10e6;
e_tol=1.0e-3;%tolerant error
phi=30/180*pi;%frictional angle
eta=6*sin(phi)/sqrt(3)/(3+sin(phi));
xi=6*cos(phi)/sqrt(3)/(3+sin(phi));
psi=0/180*pi;%dilantancy angle
eta_bar=6*sin(psi)/sqrt(3)/(3+sin(psi));

%strain increment 
inc_e=0.0001;
n_max=80;
wz=-0.5*inc_e;
%strain path
path=3;
if path==1        % oedometer test
    delta_eps=transpose([-inc_e 0 0 0 0 0 0]);
elseif path==2    % undrained triaxial compression
    delta_eps=transpose([-inc_e 0.5*inc_e 0.5*inc_e 0 0 0 0]);
elseif path==3    % simple shear
    delta_eps=transpose([0 0 0 -wz inc_e+wz 0 0]);
end

strain=zeros(7,n_max);
strain_p=zeros(n_max,1);
stress=zeros(7,n_max);
sqrt_J2=zeros(1,n_max);
stress_c=-100e3;
stress_0=[stress_c; stress_c; stress_c; 0; 0; 0; 0];
stress(:,1)=stress_0;
global stress_true strain_p_trial stress_trial;

%loop
for n=1:n_max
    delta_lamda_0=0;
    strain(:,n+1)=strain(:,n)+delta_eps;
    strain_p_trial=strain_p(n);
    stress_trial=stress(:,n)+De*delta_eps;
    phi_wave=get_residual(delta_lamda_0);
    if(phi_wave<=0)
        strain_p(n+1)=strain_p_trial;
        stress(:,n+1)=stress_trial;
        sqrt_J2(n+1)=sqrt(0.5*transpose(stress_trial)*P*stress_trial);
    else
        delta_gamma=0.1*sqrt(2/3*transpose(delta_eps)*Q*delta_eps);
        a=delta_lamda_0;
        f_a=phi_wave;
        b=delta_gamma;
        f_b=get_residual(b);
        count=0;
        while f_b>0
            b=b+delta_gamma;
            f_b=get_residual(b);
            count=count+1;
        end
        while 1
            x=a-(b-a)*f_a/(f_b-f_a);
            f_x=get_residual(x);
            if (abs(f_x)<e_tol)
                strain_p(n+1)=strain_p(n)+x;
                stress(:,n+1)=stress_true;
                sqrt_J2(n+1)=sqrt(0.5*transpose(stress_true)*P*stress_true);
                break;
            else
            
                if(f_x<0)
                    b=x;
                    f_b=f_x;
                else
                    a=x;
                    f_a=f_x;
                end
            end
        end
    end
end
p=transpose(m)*stress;

%post-process
if path==1         % oedometer test
    data.fignum=4;
    data.x=[strain(1,:);strain(3,:);strain(1,:);p];
    data.y=[stress(1,:);stress(3,:);sqrt_J2;sqrt_J2];
    data.xlabel=["$\epsilon_{xx}$";"$\epsilon_{zz}$";"$\epsilon_{xx}$";"$p$"];
    data.ylabel=["$\sigma_{xx}$/Pa";"$\sigma_{zz}$/Pa";"$\sqrt{J_2}$";"$\sqrt{J_2}$"];
    data.title=["Stress train relation";"Stress train relation";"$\epsilon_{xx}-\sqrt{J_2}$ relation";"Loading path"];
    my_plot(data)
elseif path==2
    data.fignum=4;
    data.x=[strain(1,:);strain(3,:);strain(1,:);p];
    data.y=[stress(1,:);stress(3,:);sqrt_J2;sqrt_J2];
    data.xlabel=["$\epsilon_{xx}$";"$\epsilon_{zz}$";"$\epsilon_{xx}$";"$p$"];
    data.ylabel=["$\sigma_{xx}$/Pa";"$\sigma_{zz}$/Pa";"$\sqrt{J_2}$";"$\sqrt{J_2}$"];
    data.title=["Stress train relation";"Stress train relation";"$\epsilon_{xx}-\sqrt{J_2}$ relation";"Loading path"];
    my_plot(data)
elseif path==3
    data.fignum=4;
    data.x=[strain(5,:);strain(5,:);strain(5,:);p];
    data.y=[stress(4,:);stress(5,:);sqrt_J2;sqrt_J2];
    data.xlabel=["$\epsilon_{yx}$";"$\epsilon_{yx}$";"$\epsilon_{yx}$";"$p$"];
    data.ylabel=["$\sigma_{xy}$/Pa";"$\sigma_{yx}$/Pa";"$\sqrt{J_2}$";"$\sqrt{J_2}$"];
    data.title=["Stress train relation";"Stress train relation";"$\epsilon_{yx}-\sqrt{J_2}$ relation";"Loading path"];
    my_plot(data)
end


function phi_wave=get_residual(delta_lamda)
    global P eta m xi De eta_bar K stress_trial stress_true;
    A=eye(7)+0.5/(xi*get_c(delta_lamda)-eta*transpose(m)*stress_trial+delta_lamda*eta*eta_bar*K)*delta_lamda*De*P;
    stress_true=A\(stress_trial-delta_lamda*eta_bar*De*m);
    phi_wave=sqrt(0.5*transpose(stress_true)*P*stress_true)+eta*transpose(stress_true)*m-xi*get_c(delta_lamda);
end

function c=get_c(delta_lamda)
    global H c0 strain_p_trial xi;
    c=c0+H*(strain_p_trial+xi*delta_lamda);
end