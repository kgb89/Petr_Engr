
clear,close all


mu=1 ;%cp

h=50;
nx=20; ny=20;
ts=50;

xx = [0.1 1 8 20 20 50 50 10 5 5 20 50 100 100 100 100 80 80 51 49.9];
dx = repmat(xx,20,1);

yy = [0.1;0.2;0.5;1;2;5;15;30;30;50;50;50;100;100;100;100;100;67;50;49.2];
dy = repmat(yy,1,20);

phi = ones(20,20)*0.1; phi(1,1:9) = 0.5;
kx = ones(20,20); kx(1,1:9) = 100000; 
ky = ones(20,20); ky(1,1:9) = 100000;

h = 50;
ct = 20*10^-6;
NX = length(xx);
NY = length(yy);

A = zeros(NX*NY,NX*NY);
B = zeros(NX*NY,NX*NY);

g=ones(400,1);

Vb = dx.*dy.*h;
Ax=dx.*h;
Ay=dy.*h;

dt = 1/24;


mu = ones(20,20);

fs = zeros(20);
fn=zeros(20);
fe=zeros(20);
fw=zeros(20);
 
b=ones(20,20)*5000;
br=zeros(20,20);
u=zeros(400,ts);

n=zeros(20,20);
s=zeros(20,20);
w=zeros(20,20);
e=zeros(20,20);
c=zeros(20,20);
N=zeros(380,1);
S=zeros(380,1);
W=zeros(399,1);
E=zeros(399,1);
A=zeros(400);

flag=1;
iter=0;
time=zeros(1,ts);
cumtime=zeros(1,ts);
itercount=zeros(1,ts);

for i = 1:ts
    time(:,i) = dt;
    cumtime(:,i) = sum(time(:,:));%days
    dt=dt*1.2;
end

%Iterate
for t = 1:ts;
    
    flag=1;
    iter=0;
while flag==1;
    
   
    flag=0; 
    iter=iter+1;
    itercount(:,t)=iter;
 bold(:,:)=b(:,:);
 
R=phi.*Vb.*ct/(5.615*time(:,t));

n = vec2mat(n,20);
s = vec2mat(s,20);
e = vec2mat(e,20);
w = vec2mat(w,20);
c = vec2mat(c,20);

%mu matrix
for i=1:20
    for j =1:20
mu(i,j)=1-.0001*(5000-b(i,j));
    end
end

%fs mobility coeffs kave/uave
for i=1:19
    for j =1:20
fs(i,j)=(dy(i,j)+dy(i+1,j))/(dy(i,j)/ky(i,j)+dy(i+1,j)/ky(i+1,j))/((mu(i,j)+mu(i+1,j))/2);
    end
end


i=20;
    for j =1:20
fs(i,j)=ky(i,j)/mu(i,j);
    end

%fn coeffs
for i=2:20
    for j =1:20
fn(i,j)=(dy(i,j)+dy(i-1,j))/(dy(i,j)/ky(i,j)+dy(i-1,j)/ky(i-1,j))/((mu(i,j)+mu(i-1,j))/2);
    end
end

i=1;
    for j =1:20
fn(i,j)=ky(i,j)/mu(i,j);
    end


%fe coeffs
for i=1:20
    for j =1:19
fe(i,j)=(dx(i,j)+dx(i,j+1))/(dx(i,j)/kx(i,j)+dx(i,j+1)/kx(i,j+1))/((mu(i,j)+mu(i,j+1))/2);
    end
end

j=20;
    for i =1:20
fe(i,j)=kx(i,j)/mu(i,j);
    end

%fw coeffs
for i=1:20
    for j =2:20
fw(i,j)=(dx(i,j)+dx(i,j-1))/(dx(i,j)/kx(i,j)+dx(i,j-1)/kx(i,j-1))/((mu(i,j)+mu(i,j-1))/2);
    end
end

j=1;
    for i =1:20
fw(i,j)=kx(i,j)/mu(i,j);
    end

%north coeffs
for i=2:20
    for j =1:20
n(i,j)=1.127*10^-3*fn(i,j)*Ax(i,j)/((dy(i,j)+dy(i-1,j))/2);
    end
end
n(1,:)=0;

%south coeffs
for i=1:19
    for j =1:20
s(i,j)=1.127*10^-3*fs(i,j)*Ax(i,j)/((dy(i,j)+dy(i+1,j))/2);
    end
end
s(20,:)=0;

%east coeffs
for i=1:20
    for j =1:19
e(i,j)=1.127*10^-3*fe(i,j)*Ay(i,j)/((dx(i,j)+dx(i,j+1))/2);
    end
end
e(:,20)=0;

%west coeffs
for i=1:20
    for j =2:20

w(i,j)=1.127*10^-3*fw(i,j)*Ay(i,j)/((dx(i,j)+dx(i,j-1))/2);
    end
end
w(:,1)=0;



%Center coeffs
for i=1:20
    for j =1:20
c(i,j)=-(n(i,j)+s(i,j)+w(i,j)+e(i,j)+R(i,j));
    end
end


%Create Pentadiagonal
n=transpose(n);
s=transpose(s);
w=transpose(w);
e=transpose(e);
c=transpose(c);

n=reshape(n,400,1);
s=reshape(s,400,1);
w=reshape(w,400,1);
e=reshape(e,400,1);
C=reshape(c,400,1);

for i = 21:400
N(i-20)=n(i);
end

for i = 1:380
S(i)=s(i);
end

for i = 2:400
W(i-1)=w(i);
end

for i = 1:399
E(i)=e(i);
end

C(1,1)=C(1,1)-W(1,1);
A=diag(C,0)+diag(N,-20)+diag(S,20)+diag(W,-1)+diag(E,1);
A=sparse(A);
Q(:,:)=-R.*bold;
Q(1,1) =Q(1,1)-W(1,1)*4000;
%Solver

Q=transpose(Q);
Q=reshape(Q,400,1);
Q=sparse(Q);
b=A\Q;

u(:,t) = b;

b = vec2mat(b,20);
Q = vec2mat(Q,20);

if iter==1
        flag=1;
        else
        tolerance=abs(b(:,:)-bold(:,:));
        tol=find(tolerance>1);                                      % Iteration tolerance (psi)
            if ~isempty(tol);
                flag=1;
            end
        end
end

end

% Plot Variable Grid
XGV = zeros(1,NX);
XGV(1,1) = dx(1,1);
YGV = zeros(NY,1);
YGV(1,1) = dy(1,1);
for j = 2:NX;
    XGV(1,j) = dx(1,j)+ XGV(1,j-1);
end
for j = 2:NY;
    YGV(j) = dy(j)+YGV(j-1);
end
[X,Y]=meshgrid(XGV,YGV);


t=zeros(20,20);

t0 = vec2mat(u(:,1),20);

t(:,:,1) = vec2mat(u(:,10),20);%1day
t(:,:,2) = vec2mat(u(:,18),20);%5.3
t(:,:,3) = vec2mat(u(:,22),20);%11.3
t(:,:,4) = vec2mat(u(:,25),20);%19.7
t(:,:,5) = vec2mat(u(:,27),20);%28.4
t(:,:,6) = vec2mat(u(:,29),20);%41
t(:,:,7) = vec2mat(u(:,30),20);%49.2
t(:,:,8) = vec2mat(u(:,35),20);%123
t(:,:,9) = vec2mat(u(:,14),20);%2.5

surface(X,Y,t0);
for i = 1:9
figure;
surface(X,Y,t(:,:,i));
title('Part C');
xlabel('X (ft)');
ylabel('Y (ft)');
l = colorbar;
title(l,'PSI');
end


t1 = t(:,:,1);
t2 = t(:,:,2);
t3 = t(:,:,3);
t4 = t(:,:,4);
t5 = t(:,:,5);
t6 = t(:,:,6);
t7 = t(:,:,7);
t8 = t(:,:,8);
t9 = t(:,:,9);


