% ��������ģ����ˮ������MCC�����ƺ���_MCC��������

clc;clear all;
lambda=0.0853;
kappa=0.0188;
e0=0.68;
la=lambda/(1+e0);
ka=kappa/(1+e0);
u=0.3;
p0=196; % �㶨ֵ��
% % �ƻ�Ӧ���ȣ�% �ڽ���ģ�����õ���D_Pǿ��׼��M��phi��һһ��Ӧ����ģ�
phi=35/180*pi;
M=6*sin(phi)/(3-sin(phi));
% %% �жϼ���·����1��2��3��4��5��6,7��
lujing=3
% figure
% %% ����������������
if lujing==1 % ��3=const����������ѹ������CTC��
    a1m=(3+2*M)/(3-M)*p0
    da1=1;    da3=0;    da2=0;
    a1=p0:da1:a1m;
    a3=p0*ones(size(a1));
    a2=a3;
elseif lujing==2 % P=const��ѹ������PTC��
    a1m=p0+M*p0*2/3
    da1=0.1;    da3=-0.5*da1;    da2=-0.5*da1;
    a1=p0:da1:a1m;
    a3=(3*p0-a1)/2;
    a2=a3;
elseif lujing==3 % ��1=const����ѹ����ѹ������RTC��
    a1m=(3-M)/(3+2*M)*p0
    da1=0;    da3=-0.1;    da2=-0.1;
    a3=p0:da3:a1m;
    a2=a3;
    a1=p0*ones(size(a3));
elseif lujing==4 % ��1=const����ѹ�������졪��RTE��
    a1m=(3-2*M)/(3+M)*p0
    da1=0;    da2=0;    da3=-0.1;
    a3=p0:da3:a1m;
    a1=p0*ones(size(a3));
    a2=a1;
elseif lujing==5 % P=const�����졪��PTE��
    a1m=p0+M*p0/3
    da1=0.1;    da2=da1;    da3=-2*da1;
    a1=p0:da1:a1m;
    a2=a1;
    a3=p0*3-2*a1;
elseif lujing==6 % ��3=const�������������졪��CTE��
    a1m=(3+M)/(3-2*M)*p0
    da1=01;    da2=da1;    da3=0;
    a1=p0:da1:a1m;
    a2=a1;
    a3=p0*ones(size(a1));
elseif lujing==7 % K0=1ʱΪ����̽ᣬK0�̽ᡪ��K0��
    K0=1
    a1m=10*p0;
    da1=1;    da3=da1/K0;    da2=da3;
    a1=p0:da1:a1m;
    a3=a1/K0;
    a2=a3;
end
n=length(a1);
p=1/3*(a1+a2+a3);
q=1/sqrt(2)*((a1-a2).^2+(a2-a3).^2+(a3-a1).^2).^0.5;
dp=[0,diff(p)];
dq=[0,diff(q)];
% %%----------------������Ӧ������----------------%%%%
px=q.^2/M^2./p+p; % Ӳ�����������������еı仯��
dfdp=2*p-px;
dfdq=2*q/M^2;
dfdevp=p.*px/(la-ka);
dgdp=dfdp; % �������������
dgdq=dfdq;
dlambda=(dfdp.*dp+dfdq.*dq)./dfdevp./dgdp;
dgda1=dgdp/3+dgdq*3.*(a1-p)/2./q;
dgda2=dgdp/3+dgdq*3.*(a2-p)/2./q;
dgda3=dgdp/3+dgdq*3.*(a3-p)/2./q;
dep1=dlambda.*dgda1;
dep2=dlambda.*dgda2;
dep3=dlambda.*dgda3;
% %% %%----------------����Ӧ������----------------%%%%
E=3*(1-2*u)*p/ka;                      % ���㵯��ģ��������
dee1=(1+u)./E*da1-u./E*(da1+da2+da3);
dee2=(1+u)./E*da2-u./E*(da1+da2+da3);
dee3=(1+u)./E*da3-u./E*(da1+da2+da3);
% %% �жϼ�ж�أ�
px=q.^2/M^2./p+p;
for i=1:n
    if px(i)>196
        de1(i)=dee1(i)+dep1(i);
        de2(i)=dee2(i)+dep2(i);
        de3(i)=dee3(i)+dep3(i);
    else
        de1(i)=dee1(i);
        de2(i)=dee2(i);
        de3(i)=dee3(i);
    end
end
% %% ������Ӧ�䣻
e1(1)=0;e2(1)=0;e3(1)=0;
for i=2:n
    e1(i)=e1(i-1)+de1(i);
    e2(i)=e2(i-1)+de2(i);
    e3(i)=e3(i-1)+de3(i);
end
ev=e1+e2+e3;

% %% ��ͼӦ��Ӧ���ϵ��
% figure
x1=e1*100;x2=e2*100;x3=e3*100;yev=-ev*10; eta=(q./p); % eta=(a1./a3-1); % 
plot(x1,eta,x2,eta,'-',x3,eta,x1,yev,'color','k','linewidth',2);grid on;hold on;
xlabel('��_3*100�ͦ�_1*100','fontsize',15);ylabel('��_v*10�ͦ�_1/��_3-1','fontsize',15);

