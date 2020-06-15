clc;
clear all

xx = [0.1 1 8 20 20 50 50 10 5 5 20 50 100 100 100 100 80 80 51 49.9];
dx = repmat(xx,20,1);
yy = [0.1;0.2;0.5;1;2;5;15;30;30;50;50;50;100;100;100;100;100;67;50;49.2];
dy = repmat(yy,1,20);

%porosity and permeability matrices
phi = ones(20,20)*0.1;
phi(1,1:9) = 0.5;
kx = ones(20,20);
kx(1,1:9) = 100000;
ky = ones(20,20);
ky(1,1:9) = 100000;


dt = 1/(24*5);
h = 50;
ct = 20*10^-6;
Pbc = 4000;
NX = length(xx);
NY = length(yy);

A = zeros(NX*NY,NX*NY);
R = dx.*dy.*phi.*ct/dt;

% Iterate
nstep = 7;
n = 1;
mu = ones(NX,NY,1);
mu_itr = zeros(NY,NX);       % Iteration matrix for viscosities
b = ones(NY,NX,n)*5000;
b_itr = zeros(nstep,n);      % Vector to save # iterations at each timestep
b_item = zeros(NY,NX);       % Matrix to store all iterations
b_old = ones(NY,NX)*5000;
b_new = ones(NY,NX)*5000;
b_initial = ones(NY,NX)*5000;
Rb = zeros(NY,NX);
% b_old = R.*ones(NY,NX)*5000;


for n=2:nstep+1
     dt = dt*exp(n-2);   % Variable dt
    flag = 1;
    v = 0;
    while flag == 1;
        flag = 0;
        v = v+1;
        b_old(:,:) = b_new(:,:);
        Rb(:,:) = -1*R.*b_old;
        b_item = [b_item;b_old];
        mu(:,:,v) = mu(:,:,1)-0.0001*(b_initial(:,:)-b_old(:,:));


for j = 1:NY;
    if j == 1;
        for i = 1:NX;
            if i == 1;
% harmonic permeability average, arithmetic mu average
fE = ((dx(j,i)+dx(j,i+1))/(dx(j,i)/kx(j,i)+dx(j,i+1)/kx(j,i+1)))/((mu(j,i,v)+mu(j,i+1,v))/2);
fW = kx(j,i)/mu(j,i,v);
fS = ((dy(j,i)+dy(j+1,i))/(dy(j,i)/ky(j,i)+dy(j+1,i)/ky(j+1,i)))/((mu(j,i,v)+mu(j+1,i,v))/2);
fN = ky(j,i)/mu(j,i,v);
                E = -(2*fE*(dy(j,i)/(dx(j,i+1)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j+1,i)+dy(j,i))));
                C = (E+S+W-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mE = (j-1)*NX+(i+1);
                A(mr,mE) = -E;
                mS = (j)*NX+i;
                A(mr,mS) = -S;
                Rb(1,1) = -1*R(1,1)*b_old(1,1)+Pbc*W;
            elseif i == 20;
fE = kx(j,i)/mu(j,i,v);
fW = ((dx(j,i)+dx(j,i-1))/(dx(j,i)/kx(j,i)+dx(j,i-1)/kx(j,i-1)))/((mu(j,i,v)+mu(j,i-1,v))/2);
fS = ((dy(j,i)+dy(j+1,i))/(dy(j,i)/ky(j,i)+dy(j+1,i)/ky(j+1,i)))/((mu(j,i,v)+mu(j+1,i,v))/2);
fN = ky(j,i)/mu(j,i,v);
                E = -(2*fE*(dy(j,i)/(dx(j,i)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i-1)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j+1,i)+dy(j,i))));
                C = (W+S-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mW = (j-1)*NX+(i-1);
                A(mr,mW) = -W;
                mS = (j)*NX+i;
                A(mr,mS) = -S;
            else
fE = ((dx(j,i)+dx(j,i+1))/(dx(j,i)/kx(j,i)+dx(j,i+1)/kx(j,i+1)))/((mu(j,i,v)+mu(j,i+1,v))/2);
fW = ((dx(j,i)+dx(j,i-1))/(dx(j,i)/kx(j,i)+dx(j,i-1)/kx(j,i-1)))/((mu(j,i,v)+mu(j,i-1,v))/2);
fS = ((dy(j,i)+dy(j+1,i))/(dy(j,i)/ky(j,i)+dy(j+1,i)/ky(j+1,i)))/((mu(j,i,v)+mu(j+1,i,v))/2);
fN = ky(j,i)/mu(j,i,v);
                E = -(2*fE*(dy(j,i)/(dx(j,i+1)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i-1)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j+1,i)+dy(j,i))));
                C = (E+W+S-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mE = (j-1)*NX+(i+1);
                A(mr,mE) = -E;
                mW = (j-1)*NX+(i-1);
                A(mr,mW) = -W;
                mS = (j)*NX+i;
                A(mr,mS) = -S;
            end
        end
        
    elseif j == 20;
        for i = 1:NX;
            if i == 1;
fE = ((dx(j,i)+dx(j,i+1))/(dx(j,i)/kx(j,i)+dx(j,i+1)/kx(j,i+1)))/((mu(j,i,v)+mu(j,i+1,v))/2);
fW = kx(j,i)/mu(j,i,v);
fS = ky(j,i)/mu(j,i,v);
fN = ((dy(j,i)+dy(j-1,i))/(dy(j,i)/ky(j,i)+dy(j-1,i)/ky(j-1,i)))/((mu(j,i,v)+mu(j-1,i,v))/2);
                E = -(2*fE*(dy(j,i)/(dx(j,i+1)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j-1,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j,i)+dy(j,i))));
                C = (E+N-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mE = (j-1)*NX+(i+1);
                A(mr,mE) = -E;
                mN = (j-2)*NX+i;
                A(mr,mN) = -N;
            elseif i == 20;
fE = kx(j,i)/mu(j,i,v);
fW = ((dx(j,i)+dx(j,i-1))/(dx(j,i)/kx(j,i)+dx(j,i-1)/kx(j,i-1)))/((mu(j,i,v)+mu(j,i-1,v))/2);
fS = ky(j,i)/mu(j,i,v);
fN = ((dy(j,i)+dy(j-1,i))/(dy(j,i)/ky(j,i)+dy(j-1,i)/ky(j-1,i)))/((mu(j,i,v)+mu(j-1,i,v))/2);
                E = -(2*fE*(dy(j,i)/(dx(j,i)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i-1)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j-1,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j,i)+dy(j,i))));
                C = (W+N-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mW = (j-1)*NX+(i-1);
                A(mr,mW) = -W;
                mN = (j-2)*NX+i;
                A(mr,mN) = -N;
            else
fE = ((dx(j,i)+dx(j,i+1))/(dx(j,i)/kx(j,i)+dx(j,i+1)/kx(j,i+1)))/((mu(j,i,v)+mu(j,i+1,v))/2);
fW = ((dx(j,i)+dx(j,i-1))/(dx(j,i)/kx(j,i)+dx(j,i-1)/kx(j,i-1)))/((mu(j,i,v)+mu(j,i-1,v))/2);
fS = ky(j,i)/mu(j,i,v);
fN = ((dy(j,i)+dy(j-1,i))/(dy(j,i)/ky(j,i)+dy(j-1,i)/ky(j-1,i)))/((mu(j,i,v)+mu(j-1,i,v))/2);
                E = -(2*fE*(dy(j,i)/(dx(j,i+1)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i-1)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j-1,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j,i)+dy(j,i))));
                C = (E+W+N-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mE = (j-1)*NX+(i+1);
                A(mr,mE) = -E;
                mW = (j-1)*NX+(i-1);
                A(mr,mW) = -W;
                mN = (j-2)*NX+i;
                A(mr,mN) = -N;
            end
        end
    else
        for i = 1:NX;
            if i == 1;
fE = ((dx(j,i)+dx(j,i+1))/(dx(j,i)/kx(j,i)+dx(j,i+1)/kx(j,i+1)))/((mu(j,i,v)+mu(j,i+1,v))/2);
fW = kx(j,i)/mu(j,i,v);
fS = ((dy(j,i)+dy(j+1,i))/(dy(j,i)/ky(j,i)+dy(j+1,i)/ky(j+1,i)))/((mu(j,i,v)+mu(j+1,i,v))/2);
fN = ((dy(j,i)+dy(j-1,i))/(dy(j,i)/ky(j,i)+dy(j-1,i)/ky(j-1,i)))/((mu(j,i,v)+mu(j-1,i,v))/2);
                E = -(2*fE*(dy(j,i)/(dx(j,i+1)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j-1,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j+1,i)+dy(j,i))));
                C = (E+N+S-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mE = (j-1)*NX+(i+1);
                A(mr,mE) = -E;
                mS = (j)*NX+i;
                A(mr,mS) = -S;
                mN = (j-2)*NX+i;
                A(mr,mN) = -N;
            elseif i == 20;
fE = kx(j,i)/mu(j,i,v);
fW = ((dx(j,i)+dx(j,i-1))/(dx(j,i)/kx(j,i)+dx(j,i-1)/kx(j,i-1)))/((mu(j,i,v)+mu(j,i-1,v))/2);
fS = ((dy(j,i)+dy(j+1,i))/(dy(j,i)/ky(j,i)+dy(j+1,i)/ky(j+1,i)))/((mu(j,i,v)+mu(j+1,i,v))/2);
fN = ((dy(j,i)+dy(j-1,i))/(dy(j,i)/ky(j,i)+dy(j-1,i)/ky(j-1,i)))/((mu(j,i,v)+mu(j-1,i,v))/2);
                E = -(2*fE*(dy(j,i)/(dx(j,i)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i-1)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j-1,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j+1,i)+dy(j,i))));
                C = (W+N+S-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mW = (j-1)*NX+(i-1);
                A(mr,mW) = -W;
                mS = (j)*NX+i;
                A(mr,mS) = -S;
                mN = (j-2)*NX+i;
                A(mr,mN) = -N;
            else
fE = ((dx(j,i)+dx(j,i+1))/(dx(j,i)/kx(j,i)+dx(j,i+1)/kx(j,i+1)))/((mu(j,i,v)+mu(j,i+1,v))/2);
fW = ((dx(j,i)+dx(j,i-1))/(dx(j,i)/kx(j,i)+dx(j,i-1)/kx(j,i-1)))/((mu(j,i,v)+mu(j,i-1,v))/2);
fS = ((dy(j,i)+dy(j+1,i))/(dy(j,i)/ky(j,i)+dy(j+1,i)/ky(j+1,i)))/((mu(j,i,v)+mu(j+1,i,v))/2);
fN = ((dy(j,i)+dy(j-1,i))/(dy(j,i)/ky(j,i)+dy(j-1,i)/ky(j-1,i)))/((mu(j,i,v)+mu(j-1,i,v))/2);
                E = -(2*fE*(dy(j,i)/(dx(j,i+1)+dx(j,i))));
                W = -(2*fW*(dy(j,i)/(dx(j,i-1)+dx(j,i))));
                N = -(2*fN*(dx(j,i)/(dy(j-1,i)+dy(j,i))));
                S = -(2*fS*(dx(j,i)/(dy(j+1,i)+dy(j,i))));
                C = (E+W+N+S-R(j,i));
                mr = (j-1)*NY+i;
                mC = (j-1)*NX+i;
                A(mr,mC) = C;
                mE = (j-1)*NX+(i+1);
                A(mr,mE) = -E;
                mW = (j-1)*NX+(i-1);
                A(mr,mW) = -W;
                mS = (j)*NX+i;
                A(mr,mS) = -S;
                mN = (j-2)*NX+i;
                A(mr,mN) = -N;
            end
        end
    end
end

% [r,q,s] = find(A);
% [m,g] = size(A);
% A = sparse(r,q,s,m,g);
Rb_transpose = transpose(Rb);
Rb_vec = reshape(Rb_transpose,400,1);
b_new = A\Rb_vec;
b_new = vec2mat(b_new,20);

if v == 1;              % ensures that more than one iteration is run
    flag = 1;
elseif v == 1000;
    flag = 0;
else
    tolerance = abs(b_new(:,:)-b_old(:,:));
    tol = find(tolerance>1);
    if ~isempty(tol);
        flag = 1;
    end
end
    end

b(:,:,n) = b_new(:,:);
mu(:,:,n) = mu_itr(:,:);
b_itr(n-1) = v-1;

% Plot Variable Grid
XGV = zeros(1,NX);
XGV(1,1) = xx(1,1);
YGV = zeros(NY,1);
YGV(1,1) = dy(1,1);
for i = 2:NX;
    XGV(i) = xx(i)+ XGV(i-1);
end
for j = 2:NY;
    YGV(j) = yy(j)+YGV(j-1);
end
[X,Y]=meshgrid(XGV,YGV);
Z = b(:,:,n);
caxis([4000 4220]);
title('Pressure Drop With Time');
xlabel('X (ft)');
ylabel('Y (ft)');
zlabel('Pressure (psi)');
%l = colorbar;
%title(l,'PSI');
surface(X,Y,Z)
%figure

end


