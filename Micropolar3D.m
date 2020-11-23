P1=[2 -1 -1;-1 2 -1;-1 -1 2];
P2=blkdiag([3/2 3/2; 3/2 3/2],[3/2 3/2; 3/2 3/2],[3/2 3/2; 3/2 3/2]);
P3=3*eye(9);
P=blkdiag(P1,P2,P3);

M1=eye(3);
M2=P2/3;
M3=eye(9);
M=blkdiag(M1,M2,M3);
m=zeros(18,1);
m(1)=1;m(2)=1;m(3)=1;
%m=transpose(m);

syms sg11 sg12 sg13 sg21 sg22 sg23 sg31 sg32 sg33 cp11 cp12 cp13 cp21 cp22 cp23 cp31 cp32 cp33;
syms lamda G Gc
sigma=transpose([sg11 sg22 sg33 sg12 sg21 sg23 sg32 sg31 sg13 cp11 cp22 cp33 cp12 cp13 cp21 cp23 cp31 cp32]);
p=(sg11+sg22+sg33)/3;
sigma_m=transpose([p p p 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
s=sigma-sigma_m;
Du1=lamda*ones(3,3)+2*G*eye(3);
Du2=[G+Gc G-Gc;G-Gc G+Gc];
Du=blkdiag(Du1,Du2,Du2,Du2);
Dw=2*G*eye(9);
D=blkdiag(Du,Dw);
