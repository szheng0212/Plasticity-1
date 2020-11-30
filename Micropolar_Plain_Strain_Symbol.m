a=[0.25,0.25,0.5];

P1=[2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3];
P2=[2*a(1) 2*a(2); 2*a(2) 2*a(1)];
P3=2*a(3)*eye(2);
P=blkdiag(P1,P2,P3);

M1=eye(3);
M2=P2/3;
M3=eye(9);
M=blkdiag(M1,M2,M3);
m=zeros(7,1);
m(1)=1/3;m(2)=1/3;m(3)=1/3;
%m=transpose(m);

syms sg11 sg12 sg21 sg22 sg33 cp13 cp23 eta eta_bar;
syms lamda G Gc
sigma=transpose([sg11 sg22 sg33 sg12 sg21 cp13 cp23]);
sigma_s=transpose([sg11 sg22 sg33 (sg12+sg21)/2 (sg21+sg12)/2 cp13 cp23]);
sigma_k=transpose([0 0 0 (sg12-sg21)/2 (sg21-sg12)/2 0 0]);
sg_m=(sg11+sg22+sg33)/3;
sigma_m=transpose([sg_m sg_m sg_m 0 0 0 0]);
s=sigma-sigma_m;
De1=lamda*ones(3,3)+2*G*eye(3);
De2=[G+Gc G-Gc;G-Gc G+Gc];
De3=2*G*eye(2);
De=blkdiag(De1,De2,De3);
dfdsigma=0.5*P*sigma/sqrt(0.5*transpose(sigma)*P*sigma)+ eta*m;
dgdsigma=0.5*P*sigma/sqrt(0.5*transpose(sigma)*P*sigma)+ eta_bar*m;

