%% define constants %%
clearvars

muB=9.27400968e-24; %J/T
h=6.62606957e-34; %J*s
hbar=1.054571726e-34; %J*s
e=1.602176565e-19; %C
e0=8.854187817e-12; %dielectric constant

% Energy unit is Hz (E/h)
A=117e6; %contact hyperfine
B0=0.2*2;
L1=7.5e-9; %half interdot distance
L2=7.5e-9; %half interdot distance
dg1=-0.002;%relative change in g-factor (delta_g/g)
dg2=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0];
sigma_z=[1 0; 0 -1];

% define Hamiltonian matrices %%
% basis => orbital_1 x electron_1 x nucleus_1 x orbital_2 x electron_2 x nucleus_2  

H_Hyper1=(tensor(eye(2),sigma_x,sigma_x,eye(8))+tensor(eye(2),sigma_y,sigma_y,eye(8))+tensor(eye(2),sigma_z,sigma_z,eye(8)))*A/4*tensor(eye(2)/2-sigma_x/2,eye(32));
H_Hyper2=(tensor(eye(16),sigma_x,sigma_x)+tensor(eye(16),sigma_y,sigma_y)+tensor(eye(16),sigma_z,sigma_z))*A/4*tensor(eye(8),eye(2)/2-sigma_x/2,eye(4));
H_Znuc1=-tensor(eye(4),17.2e6*B0/2*sigma_z,eye(8));
H_Znuc2=-tensor(eye(32),17.2e6*B0/2*sigma_z);
H_Zeeman1=tensor(eye(2)+(eye(2)/2+sigma_x/2)*dg1,28e9*B0/2*sigma_z,eye(16));
H_Zeeman2=tensor(eye(8),eye(2)+(eye(2)/2+sigma_x/2)*dg2,28e9*B0/2*sigma_z,eye(2));

%% 2D Edc1=Edc2=0, Vt1=Vt2=11.56e9, z space
zdim=100;
zmin=180e-9;
zmax=1260e-9;
z=zmin:(zmax-zmin)/(zdim-1):zmax;%INTERDONOR SEPARATION

gff=zeros(zdim,1);
leak=zeros(zdim,1);

Edc=-170;
H_Edc1=tensor(e*Edc*L1/h*sigma_x,eye(32));
H_Edc2=tensor(eye(8),e*Edc*L2/h*sigma_x,eye(4));
Vt=11.56e9;
H_tunel1=tensor(Vt/2*sigma_z,eye(32));
H_tunel2=tensor(eye(8),Vt/2*sigma_z,eye(4));
[Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc*L2/h*sigma_x);
singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));

for kk=1:zdim
    H_dipol=1/(h*16*pi*e0*11.7*z(kk)^3)*(e*2*L1)*(e*2*L2)*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));
    H=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_Edc1+H_Edc2+H_dipol;
    [H_eVec,H_eVal]=eig(H);
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:64
        AA=abs(singlet*H_eVec(:,i))^2;
        BB=abs(triplet*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    gff(kk)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
    leak(kk)=1-(proj1+proj2)/2;
end
%
figure(2221);
subplot(1,2,1);
%imagesc(z,Vt,log(gff')./log(10));colorbar;
plot(z/1e9,log10(gff))
ylabel('log_{10}gff (Hz)')
xlabel('z (nm)')
subplot(1,2,2);
plot(z/1e9,leak)
ylabel('leak')
xlabel('z (nm)')

%% 2D Edc1=Edc2=-2gdd/2eL, z & Vt1=Vt2 space
%Figure 5d
zdim=200;
zmin=100e-9;
zmax=500e-9;
z=zmin:(zmax-zmin)/(zdim-1):zmax;%INTERDONOR SEPARATION
Vtdim=400;
Vtmin=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)-600e6;
Vtmax=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+600e6;
Vt=Vtmin:(Vtmax-Vtmin)/(Vtdim-1):Vtmax;
gff=zeros(zdim,Vtdim);
leak=zeros(zdim,Vtdim);

for jj=1:Vtdim
    jj
    H_tunel1=tensor(Vt(jj)/2*sigma_z,eye(32));
    H_tunel2=tensor(eye(8),Vt(jj)/2*sigma_z,eye(4));
    for kk=1:zdim %distance r
        gdd=1/(h*16*pi*e0*11.7*z(kk)^3)*(e*2*L1)*(e*2*L2); %Vdd
        Edc=-2*gdd*h/(e*2*L1); %ionization point in the presence of 2nd qubit
        H_Edc1=tensor(e*Edc*L1/h*sigma_x,eye(32));
        H_Edc2=tensor(eye(8),e*Edc*L2/h*sigma_x,eye(4));
        [Orb_eVec,Orb_eVal]=eig(Vt(jj)/2*sigma_z+e*Edc*L2/h*sigma_x);
        %2 qubit states with orb in ground state and flip flop
        singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
        triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
        H_dipol=gdd*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4))); %H_dip
        H=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_Edc1+H_Edc2+H_dipol;
        [H_eVec,H_eVal]=eig(H);
        proj1=0;
        proj2=0;
        aa=0;
        bb=0;
        for i=1:64
            AA=abs(singlet*H_eVec(:,i))^2; %overlap of singlet and EigV
            BB=abs(triplet*H_eVec(:,i))^2; %overlap of triplet and EigV
            %find max
            if AA>proj1
                proj1=AA;
                aa=i;
            end
            if BB>proj2
                proj2=BB;
                bb=i;
            end
        end
        gff(kk,jj)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2); %coupling strength between singlet and triplet
        leak(kk,jj)=1-(proj1+proj2)/2;
    end
end
%
figure(222);
subplot(1,2,1);
%imagesc(z,Vt,log(gff')./log(10));colorbar;
imagesc(z,Vt,gff');axis xy;colorbar;
title('gff (Hz)')
ylabel('Vt (Hz)')
xlabel('z (m)')
subplot(1,2,2);
imagesc(Vt,z,leak');axis xy;colorbar;

%% Origin
zz=z;
VVt=Vt';
Ggff=gff';
min(min(gff))
max(max(gff))

%% z=180e-9 - 2D Edc1=Edc2, Vt1=Vt2 space
%figure 5e

z=180e-9; %distance
gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
H_dipol=gdd*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));

Edcdim=300;
Edcmin=-5000*B0;
Edcmax=5000*B0;
Edc=-2*gdd*h/(e*2*L1)+(Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax);
Vtdim=150;
Vtmin=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)-600e6;
Vtmax=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+600e6;
Vt=Vtmin:(Vtmax-Vtmin)/(Vtdim-1):Vtmax;
gff=zeros(Edcdim,Vtdim);
leak=zeros(Edcdim,Vtdim);

for jj=1:Vtdim
    jj
    H_tunel1=tensor(Vt(jj)/2*sigma_z,eye(32));
    H_tunel2=tensor(eye(8),Vt(jj)/2*sigma_z,eye(4));
    for kk=1:Edcdim
    [Orb_eVec,Orb_eVal]=eig(Vt(jj)/2*sigma_z+e*Edc(kk)*L2/h*sigma_x);
    singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
    triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
    H_Edc1=tensor(e*Edc(kk)*L1/h*sigma_x,eye(32));
    H_Edc2=tensor(eye(8),e*Edc(kk)*L2/h*sigma_x,eye(4));
    H=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_Edc1+H_Edc2+H_dipol;
    [H_eVec,H_eVal]=eig(H);
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:64
        AA=abs(singlet*H_eVec(:,i))^2;
        BB=abs(triplet*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    gff(kk,jj)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
    leak(kk,jj)=1-(proj1+proj2)/2;
    end
end
%%
figure(2221);
subplot(1,2,1);
%imagesc(z,Vt,log(gff')./log(10));
imagesc(Edc+2*gdd*h/(e*2*L1),Vt,log(gff')/log(10));axis xy;colorbar;
title('gff (Hz)')
ylabel('Vt (Hz)')
xlabel('Edc (V/m)')
subplot(1,2,2);
imagesc(Edc+2*gdd*h/(e*2*L1),Vt,leak');axis xy;colorbar;
title('leakage')
ylabel('Vt (Hz)')
xlabel('Edc (V/m)')

%% Origin
EEdc=Edc+2*gdd*h/(e*2*L1);
VVt=Vt';
Ggff=gff';
min(min(gff))
max(max(gff))

%% find adiabatic t(E) and measure gff
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+360e6; %double-dot tunel rate
z=180e-9; %INTERDONOR SEPARATION
gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
H_dipol=gdd*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));
H_tunel1=tensor(Vt/2*sigma_z,eye(32));
H_tunel2=tensor(eye(8),Vt/2*sigma_z,eye(4));
H_0=H_Znuc1+H_Znuc2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_tunel1+H_tunel2+H_dipol;

% find adiabatic t(E)
Edc00=-20000;
EdcF=-2*gdd*h/(e*2*L1);
Edcdim=2000;
dEdcMeas=(EdcF-Edc00)/(Edcdim-1);
Edc=Edc00:dEdcMeas:EdcF;

dt=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);
K=20;

gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
dDorbdEz=e*2*L1*pi/h;
gorb=Vt*pi;

for kk=1:Edcdim
    eo=sqrt(Vt^2+(e*Edc(kk)*2*L1/h+2*gdd)^2);
    pos=-(e*Edc(kk)*2*L1/h+2*gdd)/eo;
    eff=sqrt(((28e9*(1+dg1*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
    Dso(kk)=(eo-gdd*Vt/eo-eff)*pi;
    gso(kk)=A/4/sqrt(2)*Vt/eo*2*pi;
    Dorb=(e*Edc(kk)*2*L1/h+2*gdd)*pi; %IS THAT TRUE?
    dtOrb(kk)=K*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5; %IS THAT TRUE?
end

dDsodEz=-derivative(Dso)/dEdcMeas;
dgsodEz=derivative(gso)/dEdcMeas;

for kk=1:Edcdim
    dtSO(kk)=K*(dDsodEz(kk)*dEdcMeas*gso(kk)-dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
    dt(kk)=max(abs(dtOrb(kk)),abs(dtSO(kk)));
end

figure(11237);
subplot(1,2,1);
plot(Edc,dt,Edc,dtSO,Edc,dtOrb)
legend('dt','dtSO','dtOrb')
ylabel('dt (s)')
xlabel('Edc (V/m)')

% measure gff
gff=zeros(1,Edcdim);
leak=zeros(1,Edcdim);

for kk=1:Edcdim
    H_Edc1=tensor(e*Edc(kk)*L1/h*sigma_x,eye(32));
    H_Edc2=tensor(eye(8),e*Edc(kk)*L2/h*sigma_x,eye(4));
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc(kk)*L1/h*sigma_x);
    singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
    triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
    H=H_0+H_Edc1+H_Edc2;
    [H_eVec,H_eVal]=eig(H);
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:64
        AA=abs(singlet*H_eVec(:,i))^2;
        BB=abs(triplet*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    gff(kk)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
    leak(kk)=1-(proj1+proj2)/2;
end

% figure(111);
% plotyy(Edc,gff,Edc,leak);%plotyy(Edc,log(gff)/log(10),Edc,leak);

SetupRot=sum(gff.*dt)*4*pi
setupdim=Edcdim;
gatedim=2000;
tgate=(pi/2-2*SetupRot)/4/pi/gff(setupdim)*1

%% 2-q gate time evolution, Edc1=Edc2 adiabatic sweep
%figure 6a
[Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc00*L1/h*sigma_x);
UP1=tensor(Orb_eVec(:,1)',[1 0],[0 1])';
DOWN1=tensor(Orb_eVec(:,1)',[0 1],[1 0])';
UP2=UP1;
DOWN2=DOWN1;

up1=tensor(eye(2),tensor([1 0],[0 1])'*tensor([1 0],[0 1]),eye(8));
down1=tensor(eye(2),tensor([0 1],[1 0])'*tensor([0 1],[1 0]),eye(8));
int1=tensor([1 -1]'*[1 -1]./2,eye(32));
don1=tensor([1 1]'*[1 1]./2,eye(32));
up2=tensor(eye(16),tensor([1 0],[0 1])'*tensor([1 0],[0 1]));
down2=tensor(eye(16),tensor([0 1],[1 0])'*tensor([0 1],[1 0]));
int2=tensor(eye(8),[1 -1]'*[1 -1]./2,eye(4));
don2=tensor(eye(8),[1 1]'*[1 1]./2,eye(4));
flipflop1=zeros(1,2*setupdim+1+gatedim);
flipflop2=zeros(1,2*setupdim+1+gatedim);
position1=zeros(1,2*setupdim+1+gatedim);
position2=zeros(1,2*setupdim+1+gatedim);
ge1=zeros(1,2*setupdim+1+gatedim);
ge2=zeros(1,2*setupdim+1+gatedim);
Psit=zeros(64,2*setupdim+1+gatedim);

OperStep=zeros(64,64,setupdim);
Operator=eye(64);

Psi0=kron(UP1,DOWN2);
Psit(:,1)=Psi0;

for t=1:setupdim
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc(t)*L1/h*sigma_x);
    ge1(t)=Psit(:,t)'*tensor(Orb_eVec(:,2)*Orb_eVec(:,2)',eye(32))*Psit(:,t)-Psit(:,t)'*tensor(Orb_eVec(:,1)*Orb_eVec(:,1)',eye(32))*Psit(:,t);
    ge2(t)=Psit(:,t)'*tensor(eye(8),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(4))*Psit(:,t)-Psit(:,t)'*tensor(eye(8),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(4))*Psit(:,t);
    flipflop1(t)=Psit(:,t)'*up1*Psit(:,t)-Psit(:,t)'*down1*Psit(:,t);
    flipflop2(t)=Psit(:,t)'*up2*Psit(:,t)-Psit(:,t)'*down2*Psit(:,t);
    position1(t)=Psit(:,t)'*don1*Psit(:,t)-Psit(:,t)'*int1*Psit(:,t);
    position2(t)=Psit(:,t)'*don2*Psit(:,t)-Psit(:,t)'*int2*Psit(:,t);
    H_Edc1=tensor(e*Edc(t)*L1/h*sigma_x,eye(32));
    H_Edc2=tensor(eye(8),e*Edc(t)*L2/h*sigma_x,eye(4));
    H=H_0+H_Edc1+H_Edc2;
    OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
    Operator=OperStep(:,:,t)*Operator;
    Psit(:,1+t)=Operator*Psi0;
end

OperStepp=expm(-1i*2*pi*H*(tgate/gatedim));
for t=(setupdim+1):(setupdim+gatedim)
    ge1(t)=Psit(:,t)'*tensor(Orb_eVec(:,2)*Orb_eVec(:,2)',eye(32))*Psit(:,t)-Psit(:,t)'*tensor(Orb_eVec(:,1)*Orb_eVec(:,1)',eye(32))*Psit(:,t);
    ge2(t)=Psit(:,t)'*tensor(eye(8),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(4))*Psit(:,t)-Psit(:,t)'*tensor(eye(8),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(4))*Psit(:,t);
    flipflop1(t)=Psit(:,t)'*up1*Psit(:,t)-Psit(:,t)'*down1*Psit(:,t);
    flipflop2(t)=Psit(:,t)'*up2*Psit(:,t)-Psit(:,t)'*down2*Psit(:,t);
    position1(t)=Psit(:,t)'*don1*Psit(:,t)-Psit(:,t)'*int1*Psit(:,t);
    position2(t)=Psit(:,t)'*don2*Psit(:,t)-Psit(:,t)'*int2*Psit(:,t);
    Operator=OperStepp*Operator;
    Psit(:,1+t)=Operator*Psi0;
end

for t=(setupdim+gatedim+1):(2*setupdim+gatedim)
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc(2*setupdim+gatedim-t+1)*L1/h*sigma_x);
    ge1(t)=Psit(:,t)'*tensor(Orb_eVec(:,2)*Orb_eVec(:,2)',eye(32))*Psit(:,t)-Psit(:,t)'*tensor(Orb_eVec(:,1)*Orb_eVec(:,1)',eye(32))*Psit(:,t);
    ge2(t)=Psit(:,t)'*tensor(eye(8),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(4))*Psit(:,t)-Psit(:,t)'*tensor(eye(8),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(4))*Psit(:,t);
    flipflop1(t)=Psit(:,t)'*up1*Psit(:,t)-Psit(:,t)'*down1*Psit(:,t);
    flipflop2(t)=Psit(:,t)'*up2*Psit(:,t)-Psit(:,t)'*down2*Psit(:,t);
    position1(t)=Psit(:,t)'*don1*Psit(:,t)-Psit(:,t)'*int1*Psit(:,t);
    position2(t)=Psit(:,t)'*don2*Psit(:,t)-Psit(:,t)'*int2*Psit(:,t);
    Operator=OperStep(:,:,2*setupdim+gatedim-t+1)*Operator;
    Psit(:,1+t)=Operator*Psi0;
end

t=2*setupdim+gatedim+1;
ge1(t)=Psit(:,t)'*tensor(Orb_eVec(:,2)*Orb_eVec(:,2)',eye(32))*Psit(:,t)-Psit(:,t)'*tensor(Orb_eVec(:,1)*Orb_eVec(:,1)',eye(32))*Psit(:,t);
ge2(t)=Psit(:,t)'*tensor(eye(8),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(4))*Psit(:,t)-Psit(:,t)'*tensor(eye(8),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(4))*Psit(:,t);
flipflop1(t)=Psit(:,t)'*up1*Psit(:,t)-Psit(:,t)'*down1*Psit(:,t);
flipflop2(t)=Psit(:,t)'*up2*Psit(:,t)-Psit(:,t)'*down2*Psit(:,t);
position1(t)=Psit(:,t)'*don1*Psit(:,t)-Psit(:,t)'*int1*Psit(:,t);
position2(t)=Psit(:,t)'*don2*Psit(:,t)-Psit(:,t)'*int2*Psit(:,t);

% plot time evolution
%time=[tsetup (tsetup(setupdim)+tgate/gatedim):tgate/gatedim:(tsetup(setupdim)+tgate) tsetup(setupdim)+tsetup+tgate+tsetup(2) 2*tsetup(setupdim)+tgate+tsetup(3)];
time=zeros(1,2*Edcdim+gatedim+1);
for kk=2:(Edcdim+1)
    time(kk)=time(kk-1)+dt(kk-1);
end
 
for kk=(Edcdim+2):(Edcdim+gatedim+1)
    time(kk)=time(kk-1)+tgate/gatedim;
end
 
for kk=(Edcdim+gatedim+2):(2*Edcdim+gatedim+1)
    time(kk)=time(kk-1)+dt(2*Edcdim+gatedim-kk+2);
end

%% plot time evolution
Edcmid=zeros(1,gatedim);
Edcmid(:)=Edc(setupdim);
gffmid=zeros(1,gatedim);
gffmid(:)=gff(setupdim);

FigHandle = figure(569);
set(FigHandle, 'Position', [500, 100, 1200, 900]);

subplot(2,1,1);
AX = plotyy(time*1e9,[Edc Edcmid Edc(setupdim:-1:1) Edc(1)]/1000,time*1e9,[gff gffmid gff(setupdim:-1:1) gff(1)]/1e6);
set(AX(1),'Xlim',[-.3 43],'XTick',[0 20 40],'Ylim',[-20.3 .1],'YTick',[-20 -10 0],'Position', [0.45 0.65 0.5 0.3],'FontSize',50,'XTickLabel','','ticklength',[0.02 0.1])
set(AX(2),'Xlim',[-.3 43],'XTick',[0 20 40],'Ylim',[-.1 4.9],'YTick',[0 2 4],'Position', [0.35 0.65 0.5 0.3],'FontSize',50,'XTickLabel','','ticklength',[0.02 0.1])
% set(AX(2),'yscale','log')
% set(AX(2),'ytickmode','auto')
% AX(2).XTickMode = 'auto'
% AX(2).XAxisLocation = 'top'
% AX(1).Box = 'off'

subplot(2,1,2);
plot(time*1e9,flipflop1,time*1e9,flipflop2,time*1e9,position1,time*1e9,position2,time*1e9,ge1,time*1e9,ge2)
legend('flipflop1','flipflop2','position1','position2','CharEx1','CharEx2','Location','southwestoutside','boxoff')
xlabel('Time (ns)')
set(gca,'Xlim',[-.3 43],'XTick',[0 20 40],'Ylim',[-1.05 1.05],'YTick',[-1 0 1],'Position', [0.35 0.15 0.5 0.45],'FontSize',50,'ticklength',[0.02 0.1])

%% fidelities
Basis1(:,1)=UP1;
Basis1(:,2)=DOWN1;
Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
Basis2(:,1)=UP2;
Basis2(:,2)=DOWN2;
Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
Basis2q(:,1)=kron(UP1,UP2);
Basis2q(:,2)=kron(UP1,DOWN2);
Basis2q(:,3)=kron(DOWN1,UP2);
Basis2q(:,4)=kron(DOWN1,DOWN2);

phase=angle(Basis2q(:,2)'*Operator*Basis2q(:,2))-angle(Basis2q(:,4)'*Operator*Basis2q(:,4));

OperId=expm(-1i*(-phase/2*(kron(eye(2),sigma_z)+kron(sigma_z,eye(2)))))*[1 0 0 0; 0 1/sqrt(2) -1i/sqrt(2) 0; 0 -1i/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];%expm(-1i*pi/4*tensor(sigma_x,sigma_x));
OperId64=zeros(64);
for ii=1:4
    for jj=1:4
        OperId64=OperId64+OperId(ii,jj)*Basis2q(:,ii)*Basis2q(:,jj)';
    end
end

Fid16 = 0;
for ii=1:4
    for jj=1:4
        Fid16=Fid16+abs((Operator*kron(Basis1(:,ii),Basis2(:,jj)))'*(OperId64*kron(Basis1(:,ii),Basis2(:,jj))))^2/16;
    end
end
Fid16
Error = 1-Fid16

%% sweep Edc noise and average
Edcnoisedim=3;
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=100/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

FidNoise=zeros(Edcnoisedim,Edcnoisedim);
for kk=1:Edcnoisedim
    Edc1=Edc+Edcnoise(kk);
    [Orb_eVec1,Orb_eVal1]=eig(Vt/2*sigma_z+e*Edc1(1)*L1/h*sigma_x);
    UP1=tensor(Orb_eVec1(:,1)',[1 0],[0 1])';
    DOWN1=tensor(Orb_eVec1(:,1)',[0 1],[1 0])';
    for jj=1:Edcnoisedim
        tic;
        Edc2=Edc+Edcnoise(jj);
        [Orb_eVec2,Orb_eVal2]=eig(Vt/2*sigma_z+e*Edc2(1)*L2/h*sigma_x);
        UP2=tensor(Orb_eVec2(:,1)',[1 0],[0 1])';
        DOWN2=tensor(Orb_eVec2(:,1)',[0 1],[1 0])';
        
        % time evolution
        
        OperStep=zeros(64,64,setupdim);
        Operator=eye(64);
        
        for t=1:setupdim
            H_Edc1=tensor(e*Edc1(t)*L1/h*sigma_x,eye(32));
            H_Edc2=tensor(eye(8),e*Edc2(t)*L2/h*sigma_x,eye(4));
            H=H_0+H_Edc1+H_Edc2;
            OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
            Operator=OperStep(:,:,t)*Operator;
        end
        
        OperStepp=expm(-1i*2*pi*H*tgate);
        Operator=OperStepp*Operator;
        
        for t=(setupdim+gatedim+1):(2*setupdim+gatedim)
            Operator=OperStep(:,:,2*setupdim+gatedim-t+1)*Operator;
        end
        
        % fidelities
        Basis1(:,1)=UP1;
        Basis1(:,2)=DOWN1;
        Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
        Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
        Basis2(:,1)=UP2;
        Basis2(:,2)=DOWN2;
        Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
        Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
        Basis2q(:,1)=kron(UP1,UP2);
        Basis2q(:,2)=kron(UP1,DOWN2);
        Basis2q(:,3)=kron(DOWN1,UP2);
        Basis2q(:,4)=kron(DOWN1,DOWN2);

        OperId64=zeros(64);
        for oo=1:4
            for nn=1:4
                OperId64=OperId64+OperId(oo,nn)*Basis2q(:,oo)*Basis2q(:,nn)';
            end
        end
        
        FidNoise(kk,jj) = 0;
        for oo=1:4
            for nn=1:4
                FidNoise(kk,jj)=FidNoise(kk,jj)+abs((Operator*kron(Basis1(:,oo),Basis2(:,nn)))'*(OperId64*kron(Basis1(:,oo),Basis2(:,nn))))^2/16;
            end
        end
        (kk-1)/Edcnoisedim+(jj-1)/Edcnoisedim/Edcnoisedim+1/Edcnoisedim/Edcnoisedim
        toc;
    end
end
FidNoiseAv=sum(sum(FidNoise))/Edcnoisedim^2
%
figure(564)
imagesc(Edcnoise,Edcnoise,FidNoise);colorbar;

%% Sweep Vt and EdcF and calculate fidelity - Edc1=Edc2
profile on
z=180e-9; %INTERDONOR SEPARATION
H_dipol=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2)*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));%+tensor(identity,identity,identity,sigma_x,identity,identity)+tensor(identity,identity,sigma_x,identity,identity,identity)

Vtdim=17;
detmin=340e6;
detmax=420e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

EdcFdim=21;
EdcFmin=-600;
EdcFmax=400;
EdcF=(EdcFmin:(EdcFmax-EdcFmin)/(EdcFdim-1):EdcFmax);

Edc0=-20000;
Edcdim=3000; % Edc setup adiabatic sweep resolution

Edcnoisedim=5;
RMSnoise=100;
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=RMSnoise/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

ttgate = zeros(Vtdim,EdcFdim);
Fid16 = zeros(Vtdim,EdcFdim);
FidNoiseAv = zeros(Vtdim,EdcFdim);

dt=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);
K=40;

gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
dDorbdEz=e*2*L1*pi/h;

gff=zeros(1,Edcdim);

for ll = 1:Vtdim
    H_tunel1=tensor(Vt(ll)/2*sigma_z,eye(32));
    H_tunel2=tensor(eye(8),Vt(ll)/2*sigma_z,eye(4));
    H_0=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_dipol;
    gorb=Vt(ll)*pi;
    
    for mm = 1:EdcFdim
        tic;
        % find adiabatic t(E)
        dEdcMeas=(EdcF(mm)-Edc0)/(Edcdim-1);
        Edc=Edc0:dEdcMeas:EdcF(mm);
        
        for kk=1:Edcdim
            eo=sqrt(Vt(ll)^2+(e*Edc(kk)*2*L1/h+2*gdd)^2);
            pos=-(e*Edc(kk)*2*L1/h+2*gdd)/eo;
            eff=sqrt(((28e9*(1+dg1*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
            Dso(kk)=(eo-gdd*Vt(ll)/eo-eff)*pi;
            gso(kk)=A/4/sqrt(2)*Vt(ll)/eo*2*pi;
            Dorb=(e*Edc(kk)*2*L1/h+2*gdd)*pi; %IS THAT TRUE?
            dtOrb(kk)=K*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5; %IS THAT TRUE?
        end
        
        dDsodEz=-derivative(Dso)/dEdcMeas;
        dgsodEz=derivative(gso)/dEdcMeas;
        
        for kk=1:Edcdim
            dtSO(kk)=K*(dDsodEz(kk)*dEdcMeas*gso(kk)+dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
            dt(kk)=max(dtOrb(kk),dtSO(kk));
        end
        
        % measure gff
        
        for kk=1:Edcdim
            H_Edc1=tensor(e*Edc(kk)*L1/h*sigma_x,eye(32));
            H_Edc2=tensor(eye(8),e*Edc(kk)*L2/h*sigma_x,eye(4));
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc(kk)*L1/h*sigma_x);
            singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
            triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
            H=H_0+H_Edc1+H_Edc2;
            [H_eVec,H_eVal]=eig(H);
            proj1=0;
            proj2=0;
            aa=0;
            bb=0;
            for i=1:64
                AA=abs(singlet*H_eVec(:,i))^2;
                BB=abs(triplet*H_eVec(:,i))^2;
                if AA>proj1
                    proj1=AA;
                    aa=i;
                end
                if BB>proj2
                    proj2=BB;
                    bb=i;
                end
            end
            gff(kk)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
        end
        
        SetupRot=sum(gff.*dt)*4*pi;
        tgate=(pi/2-2*SetupRot)/4/pi/gff(Edcdim);
        
        if sign(tgate)>=0
            % 2-q gate time evolution
            ttgate(ll,mm) = tgate + 2*sum(dt)
        
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc0*L1/h*sigma_x);
            UP1=tensor(Orb_eVec(:,1)',[1 0],[0 1])';
            DOWN1=tensor(Orb_eVec(:,1)',[0 1],[1 0])';
            UP2=UP1;
            DOWN2=DOWN1;
            
            OperStep=zeros(64,64,Edcdim);
            Operator=eye(64);
            
            for t=1:Edcdim
                H_Edc1=tensor(e*Edc(t)*L1/h*sigma_x,eye(32));
                H_Edc2=tensor(eye(8),e*Edc(t)*L2/h*sigma_x,eye(4));
                H=H_0+H_Edc1+H_Edc2;
                OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                Operator=OperStep(:,:,t)*Operator;
            end
            
            Operator=expm(-1i*2*pi*H*tgate)*Operator;
            
            for t=(Edcdim+1):(2*Edcdim)
                Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
            end
        
            % fidelities
            Basis1(:,1)=UP1;
            Basis1(:,2)=DOWN1;
            Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
            Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
            Basis2(:,1)=UP2;
            Basis2(:,2)=DOWN2;
            Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
            Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
            Basis2q(:,1)=kron(UP1,UP2);
            Basis2q(:,2)=kron(UP1,DOWN2);
            Basis2q(:,3)=kron(DOWN1,UP2);
            Basis2q(:,4)=kron(DOWN1,DOWN2);
            
            phase=angle(Basis2q(:,2)'*Operator*Basis2q(:,2))-angle(Basis2q(:,4)'*Operator*Basis2q(:,4));
            
            OperId=expm(-1i*(-phase/2*(kron(eye(2),sigma_z)+kron(sigma_z,eye(2)))))*[1 0 0 0; 0 1/sqrt(2) -1i/sqrt(2) 0; 0 -1i/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];%expm(-1i*pi/4*tensor(sigma_x,sigma_x));
            OperId64=zeros(64);
            for ii=1:4
                for jj=1:4
                    OperId64=OperId64+OperId(ii,jj)*Basis2q(:,ii)*Basis2q(:,jj)';
                end
            end
            
            Fid16(ll,mm) = 0;
            for ii=1:4
                for jj=1:4
                    Fid16(ll,mm)=Fid16(ll,mm)+abs((Operator*kron(Basis1(:,ii),Basis2(:,jj)))'*(OperId64*kron(Basis1(:,ii),Basis2(:,jj))))^2/16;
                end
            end
            Fid16(ll,mm)=Fid16(ll,mm)
        
            % sweep Edc noise and average
            FidNoise=zeros(Edcnoisedim,Edcnoisedim);
            for kk=1:Edcnoisedim
                Edc1=Edc+Edcnoise(kk);
                [Orb_eVec1,Orb_eVal1]=eig(Vt(ll)/2*sigma_z+e*Edc1(1)*L1/h*sigma_x);
                UP1=tensor(Orb_eVec1(:,1)',[1 0],[0 1])';
                DOWN1=tensor(Orb_eVec1(:,1)',[0 1],[1 0])';
                for jj=1:Edcnoisedim
                    Edc2=Edc+Edcnoise(jj);
                    [Orb_eVec2,Orb_eVal2]=eig(Vt(ll)/2*sigma_z+e*Edc2(1)*L2/h*sigma_x);
                    UP2=tensor(Orb_eVec2(:,1)',[1 0],[0 1])';
                    DOWN2=tensor(Orb_eVec2(:,1)',[0 1],[1 0])';
                    
                    % time evolution
                    
                    OperStep=zeros(64,64,Edcdim);
                    Operator=eye(64);
                    
                    for t=1:Edcdim
                        H_Edc1=tensor(e*Edc1(t)*L1/h*sigma_x,eye(32));
                        H_Edc2=tensor(eye(8),e*Edc2(t)*L2/h*sigma_x,eye(4));
                        H=H_0+H_Edc1+H_Edc2;
                        OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                        Operator=OperStep(:,:,t)*Operator;
                    end
                    
                    Operator=expm(-1i*2*pi*H*tgate)*Operator;
                    
                    for t=(Edcdim+1):(2*Edcdim)
                        Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                    end
                    
                    % fidelities
                    Basis1(:,1)=UP1;
                    Basis1(:,2)=DOWN1;
                    Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
                    Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
                    Basis2(:,1)=UP2;
                    Basis2(:,2)=DOWN2;
                    Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
                    Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
                    Basis2q(:,1)=kron(UP1,UP2);
                    Basis2q(:,2)=kron(UP1,DOWN2);
                    Basis2q(:,3)=kron(DOWN1,UP2);
                    Basis2q(:,4)=kron(DOWN1,DOWN2);
                    
                    OperId64=zeros(64);
                    for oo=1:4
                        for nn=1:4
                            OperId64=OperId64+OperId(oo,nn)*Basis2q(:,oo)*Basis2q(:,nn)';
                        end
                    end
                    
                    FidNoise(kk,jj) = 0;
                    for oo=1:4
                        for nn=1:4
                            FidNoise(kk,jj)=FidNoise(kk,jj)+abs((Operator*kron(Basis1(:,oo),Basis2(:,nn)))'*(OperId64*kron(Basis1(:,oo),Basis2(:,nn))))^2/16;
                        end
                    end
                    
                end
            end
            FidNoiseAv(ll,mm)=sum(sum(FidNoise))/Edcnoisedim^2
        else
            ttgate(ll,mm) = NaN;
            Fid16(ll,mm) = NaN;
            FidNoiseAv(ll,mm) = NaN;
        end
        save 2QFid_noise_adiabatic_05.mat
        toc;
        (ll-1)/Vtdim+(mm-1)/Vtdim/EdcFdim+1/Vtdim/EdcFdim
    end
end
profile viewer
%%
figure(678);
subplot(1,3,1);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,ttgate*1e9);axis xy;colorbar;
title('Gate time (ns)')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
% xlim([155 620])
subplot(1,3,2);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid16))/log(10));axis xy;colorbar;
title('log_{10}(Bare gate error)')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
% xlim([155 620])
subplot(1,3,3);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv)/log(10));axis xy;colorbar;
title('log_{10}(Gate error with noise)')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
% xlim([155 620])

%% Sweep K, Vt and Edc1=Edc2 around EdcF = -2*gdd*h/(e*2*L1) and calculate best fidelity for each V_t
clearvars

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C
e0=8.854187817e-12; %dielectric constant

% Energy unit is Hz (E/h)
A=117e6; %contact hyperfine
B0=0.2*2;
L1=7.5e-9; %half interdot distance
L2=7.5e-9; %half interdot distance
dg1=-0.002;%relative change in g-factor (delta_g/g)
dg2=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0];
sigma_z=[1 0; 0 -1];

% define Hamiltonian matrices %%
% basis => orbital_1 x electron_1 x nucleus_1 x orbital_2 x electron_2 x nucleus_2

H_Hyper1=(tensor(eye(2),sigma_x,sigma_x,eye(8))+tensor(eye(2),sigma_y,sigma_y,eye(8))+tensor(eye(2),sigma_z,sigma_z,eye(8)))*A/4*tensor(eye(2)/2-sigma_x/2,eye(32));
H_Hyper2=(tensor(eye(16),sigma_x,sigma_x)+tensor(eye(16),sigma_y,sigma_y)+tensor(eye(16),sigma_z,sigma_z))*A/4*tensor(eye(8),eye(2)/2-sigma_x/2,eye(4));
H_Znuc1=-tensor(eye(4),17.2e6*B0/2*sigma_z,eye(8));
H_Znuc2=-tensor(eye(32),17.2e6*B0/2*sigma_z);
H_Zeeman1=tensor(eye(2)+(eye(2)/2+sigma_x/2)*dg1,28e9*B0/2*sigma_z,eye(16));
H_Zeeman2=tensor(eye(8),eye(2)+(eye(2)/2+sigma_x/2)*dg2,28e9*B0/2*sigma_z,eye(2));
profile on
z=180e-9; %INTERDONOR SEPARATION
H_dipol=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2)*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));%+tensor(identity,identity,identity,sigma_x,identity,identity)+tensor(identity,identity,sigma_x,identity,identity,identity)

Vtdim=1;
%detmin=375e6;
%detmax=375e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+410e6;%(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

Edc0=-20000/B0;
Edcdim=2000; % Edc setup adiabatic sweep resolution

Edcnoisedim=3;
RMSnoise=100;
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=RMSnoise/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

dt=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);
Kdim=10;
Kmin=13;
Kmax=31;
K=Kmin:(Kmax-Kmin)/(Kdim-1):Kmax;

gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
dDorbdEz=e*2*L1*pi/h;

gff=zeros(1,Edcdim);

EdcFdim=13;
EdcFmax=300;
EdcF=-2*gdd*h/(e*2*L1)+(-EdcFmax:2*EdcFmax/(EdcFdim-1):EdcFmax);

ttgate = zeros(Vtdim,EdcFdim,Kdim);
Fid16 = zeros(Vtdim,EdcFdim,Kdim);
FidNoiseAv = zeros(Vtdim,EdcFdim,Kdim);
RelErr = zeros(Vtdim,EdcFdim,Kdim);

for ll = 1:Vtdim
    tic;
    H_tunel1=tensor(Vt(ll)/2*sigma_z,eye(32));
    H_tunel2=tensor(eye(8),Vt(ll)/2*sigma_z,eye(4));
    H_0=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_dipol;
    gorb=Vt(ll)*pi;
    
    for eee=1:EdcFdim
        dEdcMeas=(EdcF(eee)-Edc0)/(Edcdim-1);
        Edc=Edc0:dEdcMeas:EdcF(eee);
        
        for kadia=1:Kdim
            % find adiabatic t(E)
            for kk=1:Edcdim
                eo=sqrt(Vt(ll)^2+(e*Edc(kk)*2*L1/h+2*gdd)^2);
                pos=-(e*Edc(kk)*2*L1/h+2*gdd)/eo;
                eff=sqrt(((28e9*(1+dg1*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
                Dso(kk)=(eo-gdd*Vt(ll)/eo-eff)*pi;
                gso(kk)=A/4/sqrt(2)*Vt(ll)/eo*2*pi;
                Dorb=(e*Edc(kk)*2*L1/h+2*gdd)*pi; %IS THAT TRUE?
                dtOrb(kk)=K(kadia)*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5; %IS THAT TRUE?
            end
            
            dDsodEz=-derivative(Dso)/dEdcMeas;
            dgsodEz=derivative(gso)/dEdcMeas;
            
            for kk=1:Edcdim
                dtSO(kk)=K(kadia)*abs(dDsodEz(kk)*dEdcMeas*gso(kk)-dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
                dt(kk)=max(dtOrb(kk),dtSO(kk));
            end
            
            charup1=zeros(64,64,Edcdim);
            charup2=zeros(64,64,Edcdim);
            % measure gff
            for kk=1:Edcdim
                H_Edc1=tensor(e*Edc(kk)*L1/h*sigma_x,eye(32));
                H_Edc2=tensor(eye(8),e*Edc(kk)*L2/h*sigma_x,eye(4));
                [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc(kk)*L1/h*sigma_x);
                charup1(:,:,kk)=tensor(Orb_eVec(:,2)*Orb_eVec(:,2)',eye(32));
                charup2(:,:,kk)=tensor(eye(8),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(4));
                singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
                triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
                H=H_0+H_Edc1+H_Edc2;
                [H_eVec,H_eVal]=eig(H);
                proj1=0;
                proj2=0;
                aa=0;
                bb=0;
                for i=1:64
                    AA=abs(singlet*H_eVec(:,i))^2;
                    BB=abs(triplet*H_eVec(:,i))^2;
                    if AA>proj1
                        proj1=AA;
                        aa=i;
                    end
                    if BB>proj2
                        proj2=BB;
                        bb=i;
                    end
                end
                gff(kk)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
            end
            
            SetupRot=sum(gff.*dt)*4*pi;
            tgate=(pi/2-2*SetupRot)/4/pi/gff(Edcdim);
            
            if sign(tgate)>=0
                % 2-q gate time evolution
                ttgate(ll,eee,kadia) = tgate + 2*sum(dt);
                'gate time'
                ttgate(ll,eee,kadia)
                
                [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc0*L1/h*sigma_x);
                UP1=tensor(Orb_eVec(:,1)',[1 0],[0 1])';
                DOWN1=tensor(Orb_eVec(:,1)',[0 1],[1 0])';
                UP2=UP1;
                DOWN2=DOWN1;
                
                OperStep=zeros(64,64,Edcdim);
                Operator=eye(64);
                ge1=zeros(64,64);
                ge2=zeros(64,64);
                for t=1:Edcdim
                    H_Edc1=tensor(e*Edc(t)*L1/h*sigma_x,eye(32));
                    H_Edc2=tensor(eye(8),e*Edc(t)*L2/h*sigma_x,eye(4));
                    H=H_0+H_Edc1+H_Edc2;
                    OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                    Operator=OperStep(:,:,t)*Operator;
                    ge1=ge1+Operator'*charup1(:,:,t)*Operator*dt(t);
                    ge2=ge2+Operator'*charup2(:,:,t)*Operator*dt(t);
                end
                
                Operator=expm(-1i*2*pi*H*tgate)*Operator;
                ge1=ge1+Operator'*charup1(:,:,Edcdim)*Operator*tgate;
                ge2=ge2+Operator'*charup2(:,:,Edcdim)*Operator*tgate;
                
                for t=(Edcdim+1):(2*Edcdim)
                    Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                    ge1=ge1+Operator'*charup1(:,:,2*Edcdim-t+1)*Operator*dt(2*Edcdim-t+1);
                    ge2=ge2+Operator'*charup2(:,:,2*Edcdim-t+1)*Operator*dt(2*Edcdim-t+1);
                end
                
                % fidelities
                Basis1(:,1)=UP1;
                Basis1(:,2)=DOWN1;
                Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
                Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
                Basis2(:,1)=UP2;
                Basis2(:,2)=DOWN2;
                Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
                Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
                Basis2q(:,1)=kron(UP1,UP2);
                Basis2q(:,2)=kron(UP1,DOWN2);
                Basis2q(:,3)=kron(DOWN1,UP2);
                Basis2q(:,4)=kron(DOWN1,DOWN2);
                
                phase=angle(Basis2q(:,2)'*Operator*Basis2q(:,2))-angle(Basis2q(:,4)'*Operator*Basis2q(:,4));
                
                OperId=expm(-1i*(-phase/2*(kron(eye(2),sigma_z)+kron(sigma_z,eye(2)))))*[1 0 0 0; 0 1/sqrt(2) -1i/sqrt(2) 0; 0 -1i/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];%expm(-1i*pi/4*tensor(sigma_x,sigma_x));
                OperId64=zeros(64);
                for ii=1:4
                    for jj=1:4
                        OperId64=OperId64+OperId(ii,jj)*Basis2q(:,ii)*Basis2q(:,jj)';
                    end
                end
                
                for ii=1:4
                    for jj=1:4
                        Fid16(ll,eee,kadia)=Fid16(ll,eee,kadia)+abs((Operator*kron(Basis1(:,ii),Basis2(:,jj)))'*(OperId64*kron(Basis1(:,ii),Basis2(:,jj))))^2/16;
                    end
                end
                'Leakage error'
                1-Fid16(ll,eee,kadia)
                
                RelErr1=0;
                RelErr2=0;
                for ii=1:4
                    for jj=1:4
                        RelErr1=RelErr1+(kron(Basis1(:,ii),Basis2(:,jj))'*ge1*kron(Basis1(:,ii),Basis2(:,jj)))/16;
                        RelErr2=RelErr2+(kron(Basis1(:,ii),Basis2(:,jj))'*ge2*kron(Basis1(:,ii),Basis2(:,jj)))/16;
                    end
                end
                RelErr1=RelErr1*2.37e-24*eo*Vt(ll)^2;
                RelErr2=RelErr2*2.37e-24*eo*Vt(ll)^2;
                RelErr(ll,eee,kadia)=real(1-exp(-RelErr1))/2+real(1-exp(-RelErr2))/2;
                'Relaxation error'
                RelErr(ll,eee,kadia)
                
                % sweep Edc noise and average
                FidNoise=zeros(Edcnoisedim,Edcnoisedim);
                for kk=1:Edcnoisedim
                    Edc1=Edc+Edcnoise(kk);
                    [Orb_eVec1,Orb_eVal1]=eig(Vt(ll)/2*sigma_z+e*Edc1(1)*L1/h*sigma_x);
                    UP1=tensor(Orb_eVec1(:,1)',[1 0],[0 1])';
                    DOWN1=tensor(Orb_eVec1(:,1)',[0 1],[1 0])';
                    for jj=1:Edcnoisedim
                        Edc2=Edc+Edcnoise(jj);
                        [Orb_eVec2,Orb_eVal2]=eig(Vt(ll)/2*sigma_z+e*Edc2(1)*L2/h*sigma_x);
                        UP2=tensor(Orb_eVec2(:,1)',[1 0],[0 1])';
                        DOWN2=tensor(Orb_eVec2(:,1)',[0 1],[1 0])';
                        
                        % time evolution
                        
                        OperStep=zeros(64,64,Edcdim);
                        Operator=eye(64);
                        
                        for t=1:Edcdim
                            H_Edc1=tensor(e*Edc1(t)*L1/h*sigma_x,eye(32));
                            H_Edc2=tensor(eye(8),e*Edc2(t)*L2/h*sigma_x,eye(4));
                            H=H_0+H_Edc1+H_Edc2;
                            OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                            Operator=OperStep(:,:,t)*Operator;
                        end
                        
                        Operator=expm(-1i*2*pi*H*tgate)*Operator;
                        
                        for t=(Edcdim+1):(2*Edcdim)
                            Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                        end
                        
                        % fidelities
                        Basis1(:,1)=UP1;
                        Basis1(:,2)=DOWN1;
                        Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
                        Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
                        Basis2(:,1)=UP2;
                        Basis2(:,2)=DOWN2;
                        Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
                        Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
                        Basis2q(:,1)=kron(UP1,UP2);
                        Basis2q(:,2)=kron(UP1,DOWN2);
                        Basis2q(:,3)=kron(DOWN1,UP2);
                        Basis2q(:,4)=kron(DOWN1,DOWN2);
                        
                        OperId64=zeros(64);
                        for oo=1:4
                            for nn=1:4
                                OperId64=OperId64+OperId(oo,nn)*Basis2q(:,oo)*Basis2q(:,nn)';
                            end
                        end
                        
                        FidNoise(kk,jj) = 0;
                        for oo=1:4
                            for nn=1:4
                                FidNoise(kk,jj)=FidNoise(kk,jj)+abs((Operator*kron(Basis1(:,oo),Basis2(:,nn)))'*(OperId64*kron(Basis1(:,oo),Basis2(:,nn))))^2/16;
                            end
                        end
                        
                    end
                end
                FidNoiseAv(ll,eee,kadia)=sum(sum(FidNoise))/Edcnoisedim^2;
                'Noise error'
                1-FidNoiseAv(ll,eee,kadia)
            else
                ttgate(ll,eee,kadia) = NaN;
                Fid16(ll,eee,kadia) = NaN;
                FidNoiseAv(ll,eee,kadia) = NaN;
                RelErr(ll,eee,kadia) = NaN;
            end
            kadia,eee,ll
        end
    end
    %save 2QFid_noise_adiabatic_Vt_optimizeEdcF&K_04.mat
    toc;
    %         (ll-1)/Vtdim+(mm-1)/Vtdim/EdcFdim+1/Vtdim/EdcFdim
end
profile viewer

%%
figure(673815);
subplot(1,5,1);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv(:,:,1))/log(10));axis xy;colorbar;
title('Noise error, K=20')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
subplot(1,5,2);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv(:,:,2))/log(10));axis xy;colorbar;
title('Noise error, K=25')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
subplot(1,5,3);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv(:,:,3))/log(10));axis xy;colorbar;
title('Noise error, K=30')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
subplot(1,5,4);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv(:,:,4))/log(10));axis xy;colorbar;
title('Noise error, K=35')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')
subplot(1,5,5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv(:,:,5))/log(10));axis xy;colorbar;
title('Noise error, K=40')
xlabel('Edc final (V/m)')
ylabel('\delta_{so} (MHz)')

%%
[FidNoiseAv1,kopt] = max(FidNoiseAv,[],3);
[FidNoiseAv2,EdcFopt] = max(FidNoiseAv1,[],2);
Kopt=kopt(EdcFopt);
for i=1:Vtdim
    Kopt(i) = kopt(i,EdcFopt(i));
end
for i=1:Vtdim
    ttgate2(i) = ttgate(i,EdcFopt(i),Kopt(i));
end
for i=1:Vtdim
    Fid162(i) = Fid16(i,EdcFopt(i),Kopt(i));
end
for i=1:Vtdim
    RelErr2(i) = RelErr(i,EdcFopt(i),Kopt(i));
end

%% Origin
VVt=Vt'/1e9;
KKopt=K(Kopt)';
EEdcFopt=-EdcF(EdcFopt)';
Gtime=ttgate2'*1e9;
LeakN=1-Fid162';
TotN=1-FidNoiseAv2;
RRelErr=RelErr2';


    %%
figure(673819);
subplot(2,3,1);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,ttgate2*1e9);
ylabel('Gate time (ns)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,3,2);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,K(Kopt));
ylabel('optimal K')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,3,3);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,-EdcF(EdcFopt));
ylabel('optimal EdcF (V/m)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,3,4);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid162))/log(10));
ylabel('log_{10}(Leakage error)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,3,5);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv2)/log(10));
ylabel('log_{10}(Noise error)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,3,6);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(RelErr2)/log(10));
ylabel('log_{10}(Relaxation error)')
xlabel('V_t - \epsilon_{ff} (MHz)')


%% Sweep Vt and calculate fidelity - Edc1=Edc2, EdcF = -2*gdd*h/(e*2*L1)
clearvars

muB=9.27400968e-24; %J/T
h=6.62606957e-34; %J*s
hbar=1.054571726e-34; %J*s
e=1.602176565e-19; %C
e0=8.854187817e-12; %dielectric constant

% Energy unit is Hz (E/h)
A=117e6; %contact hyperfine
B0=0.2*2;
L1=7.5e-9; %half interdot distance
L2=7.5e-9; %half interdot distance
dg1=-0.002;%relative change in g-factor (delta_g/g)
dg2=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0];
sigma_z=[1 0; 0 -1];

% define Hamiltonian matrices %%
% basis => orbital_1 x electron_1 x nucleus_1 x orbital_2 x electron_2 x nucleus_2  

H_Hyper1=(tensor(eye(2),sigma_x,sigma_x,eye(8))+tensor(eye(2),sigma_y,sigma_y,eye(8))+tensor(eye(2),sigma_z,sigma_z,eye(8)))*A/4*tensor(eye(2)/2-sigma_x/2,eye(32));
H_Hyper2=(tensor(eye(16),sigma_x,sigma_x)+tensor(eye(16),sigma_y,sigma_y)+tensor(eye(16),sigma_z,sigma_z))*A/4*tensor(eye(8),eye(2)/2-sigma_x/2,eye(4));
H_Znuc1=-tensor(eye(4),17.2e6*B0/2*sigma_z,eye(8));
H_Znuc2=-tensor(eye(32),17.2e6*B0/2*sigma_z);
H_Zeeman1=tensor(eye(2)+(eye(2)/2+sigma_x/2)*dg1,28e9*B0/2*sigma_z,eye(16));
H_Zeeman2=tensor(eye(8),eye(2)+(eye(2)/2+sigma_x/2)*dg2,28e9*B0/2*sigma_z,eye(2));
profile on
z=180e-9; %INTERDONOR SEPARATION
H_dipol=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2)*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));%+tensor(identity,identity,identity,sigma_x,identity,identity)+tensor(identity,identity,sigma_x,identity,identity,identity)

Vtdim=24;
detmin=365e6;
detmax=480e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

Edc0=-20000/B0;
Edcdim=2000; % Edc setup adiabatic sweep resolution

Edcnoisedim=3;
RMSnoise=100;
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=RMSnoise/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

ttgate = zeros(Vtdim,1);
Fid16 = zeros(Vtdim,1);
FidNoiseAv = zeros(Vtdim,1);
RelErr = zeros(Vtdim,1);

dt=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);
K=25;

gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
dDorbdEz=e*2*L1*pi/h;

gff=zeros(1,Edcdim);

EdcF=-2*gdd*h/(e*2*L1);
dEdcMeas=(EdcF-Edc0)/(Edcdim-1);
Edc=Edc0:dEdcMeas:EdcF;

for ll = 1:Vtdim
    tic;
    H_tunel1=tensor(Vt(ll)/2*sigma_z,eye(32));
    H_tunel2=tensor(eye(8),Vt(ll)/2*sigma_z,eye(4));
    H_0=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_dipol;
    gorb=Vt(ll)*pi;
        
        % find adiabatic t(E)
        for kk=1:Edcdim
            eo=sqrt(Vt(ll)^2+(e*Edc(kk)*2*L1/h+2*gdd)^2);
            pos=-(e*Edc(kk)*2*L1/h+2*gdd)/eo;
            eff=sqrt(((28e9*(1+dg1*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
            Dso(kk)=(eo-gdd*Vt(ll)/eo-eff)*pi;
            gso(kk)=A/4/sqrt(2)*Vt(ll)/eo*2*pi;
            Dorb=(e*Edc(kk)*2*L1/h+2*gdd)*pi; %IS THAT TRUE?
            dtOrb(kk)=K*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5; %IS THAT TRUE?
        end
        
        dDsodEz=-derivative(Dso)/dEdcMeas;
        dgsodEz=derivative(gso)/dEdcMeas;
        
        for kk=1:Edcdim
            dtSO(kk)=K*abs(dDsodEz(kk)*dEdcMeas*gso(kk)-dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
            dt(kk)=max(dtOrb(kk),dtSO(kk));
        end
        
        charup1=zeros(64,64,Edcdim);
        charup2=zeros(64,64,Edcdim);
        % measure gff
        for kk=1:Edcdim
            H_Edc1=tensor(e*Edc(kk)*L1/h*sigma_x,eye(32));
            H_Edc2=tensor(eye(8),e*Edc(kk)*L2/h*sigma_x,eye(4));
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc(kk)*L1/h*sigma_x);
            charup1(:,:,kk)=tensor(Orb_eVec(:,2)*Orb_eVec(:,2)',eye(32));
            charup2(:,:,kk)=tensor(eye(8),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(4));
            singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
            triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
            H=H_0+H_Edc1+H_Edc2;
            [H_eVec,H_eVal]=eig(H);
            proj1=0;
            proj2=0;
            aa=0;
            bb=0;
            for i=1:64
                AA=abs(singlet*H_eVec(:,i))^2;
                BB=abs(triplet*H_eVec(:,i))^2;
                if AA>proj1
                    proj1=AA;
                    aa=i;
                end
                if BB>proj2
                    proj2=BB;
                    bb=i;
                end
            end
            gff(kk)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
        end
        
        SetupRot=sum(gff.*dt)*4*pi;
        tgate=(pi/2-2*SetupRot)/4/pi/gff(Edcdim);
        
        if sign(tgate)>=0
            % 2-q gate time evolution
            ttgate(ll) = tgate + 2*sum(dt)
        
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc0*L1/h*sigma_x);
            UP1=tensor(Orb_eVec(:,1)',[1 0],[0 1])';
            DOWN1=tensor(Orb_eVec(:,1)',[0 1],[1 0])';
            UP2=UP1;
            DOWN2=DOWN1;
            
            OperStep=zeros(64,64,Edcdim);
            Operator=eye(64);
            ge1=zeros(64,64);
            ge2=zeros(64,64);
            for t=1:Edcdim
                H_Edc1=tensor(e*Edc(t)*L1/h*sigma_x,eye(32));
                H_Edc2=tensor(eye(8),e*Edc(t)*L2/h*sigma_x,eye(4));
                H=H_0+H_Edc1+H_Edc2;
                OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                Operator=OperStep(:,:,t)*Operator;
                ge1=ge1+Operator'*charup1(:,:,t)*Operator*dt(t);
                ge2=ge2+Operator'*charup2(:,:,t)*Operator*dt(t);
            end
            
            Operator=expm(-1i*2*pi*H*tgate)*Operator;
            ge1=ge1+Operator'*charup1(:,:,Edcdim)*Operator*tgate;
            ge2=ge2+Operator'*charup2(:,:,Edcdim)*Operator*tgate;
            
            for t=(Edcdim+1):(2*Edcdim)
                Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                ge1=ge1+Operator'*charup1(:,:,2*Edcdim-t+1)*Operator*dt(2*Edcdim-t+1);
                ge2=ge2+Operator'*charup2(:,:,2*Edcdim-t+1)*Operator*dt(2*Edcdim-t+1);
            end
        
            % fidelities
            Basis1(:,1)=UP1;
            Basis1(:,2)=DOWN1;
            Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
            Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
            Basis2(:,1)=UP2;
            Basis2(:,2)=DOWN2;
            Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
            Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
            Basis2q(:,1)=kron(UP1,UP2);
            Basis2q(:,2)=kron(UP1,DOWN2);
            Basis2q(:,3)=kron(DOWN1,UP2);
            Basis2q(:,4)=kron(DOWN1,DOWN2);
            
            phase=angle(Basis2q(:,2)'*Operator*Basis2q(:,2))-angle(Basis2q(:,4)'*Operator*Basis2q(:,4));
            
            OperId=expm(-1i*(-phase/2*(kron(eye(2),sigma_z)+kron(sigma_z,eye(2)))))*[1 0 0 0; 0 1/sqrt(2) -1i/sqrt(2) 0; 0 -1i/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];%expm(-1i*pi/4*tensor(sigma_x,sigma_x));
            OperId64=zeros(64);
            for ii=1:4
                for jj=1:4
                    OperId64=OperId64+OperId(ii,jj)*Basis2q(:,ii)*Basis2q(:,jj)';
                end
            end
            
            for ii=1:4
                for jj=1:4
                    Fid16(ll)=Fid16(ll)+abs((Operator*kron(Basis1(:,ii),Basis2(:,jj)))'*(OperId64*kron(Basis1(:,ii),Basis2(:,jj))))^2/16;
                end
            end
            Fid16(ll)=Fid16(ll)
            
            RelErr1=0;
            RelErr2=0;
            for ii=1:4
                for jj=1:4
                    RelErr1=RelErr1+(kron(Basis1(:,ii),Basis2(:,jj))'*ge1*kron(Basis1(:,ii),Basis2(:,jj)))/16;
                    RelErr2=RelErr2+(kron(Basis1(:,ii),Basis2(:,jj))'*ge2*kron(Basis1(:,ii),Basis2(:,jj)))/16;
                end
            end
            RelErr1=RelErr1*2.37e-24*eo*Vt(ll)^2;
            RelErr2=RelErr2*2.37e-24*eo*Vt(ll)^2;
            RelErr(ll)=real(1-exp(-RelErr1))/2+real(1-exp(-RelErr2))/2
        
            % sweep Edc noise and average
            FidNoise=zeros(Edcnoisedim,Edcnoisedim);
            for kk=1:Edcnoisedim
                Edc1=Edc+Edcnoise(kk);
                [Orb_eVec1,Orb_eVal1]=eig(Vt(ll)/2*sigma_z+e*Edc1(1)*L1/h*sigma_x);
                UP1=tensor(Orb_eVec1(:,1)',[1 0],[0 1])';
                DOWN1=tensor(Orb_eVec1(:,1)',[0 1],[1 0])';
                for jj=1:Edcnoisedim
                    Edc2=Edc+Edcnoise(jj);
                    [Orb_eVec2,Orb_eVal2]=eig(Vt(ll)/2*sigma_z+e*Edc2(1)*L2/h*sigma_x);
                    UP2=tensor(Orb_eVec2(:,1)',[1 0],[0 1])';
                    DOWN2=tensor(Orb_eVec2(:,1)',[0 1],[1 0])';
                    
                    % time evolution
                    
                    OperStep=zeros(64,64,Edcdim);
                    Operator=eye(64);
                    
                    for t=1:Edcdim
                        H_Edc1=tensor(e*Edc1(t)*L1/h*sigma_x,eye(32));
                        H_Edc2=tensor(eye(8),e*Edc2(t)*L2/h*sigma_x,eye(4));
                        H=H_0+H_Edc1+H_Edc2;
                        OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                        Operator=OperStep(:,:,t)*Operator;
                    end
                    
                    Operator=expm(-1i*2*pi*H*tgate)*Operator;
                    
                    for t=(Edcdim+1):(2*Edcdim)
                        Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                    end
                    
                    % fidelities
                    Basis1(:,1)=UP1;
                    Basis1(:,2)=DOWN1;
                    Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
                    Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
                    Basis2(:,1)=UP2;
                    Basis2(:,2)=DOWN2;
                    Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
                    Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
                    Basis2q(:,1)=kron(UP1,UP2);
                    Basis2q(:,2)=kron(UP1,DOWN2);
                    Basis2q(:,3)=kron(DOWN1,UP2);
                    Basis2q(:,4)=kron(DOWN1,DOWN2);
                    
                    OperId64=zeros(64);
                    for oo=1:4
                        for nn=1:4
                            OperId64=OperId64+OperId(oo,nn)*Basis2q(:,oo)*Basis2q(:,nn)';
                        end
                    end
                    
                    FidNoise(kk,jj) = 0;
                    for oo=1:4
                        for nn=1:4
                            FidNoise(kk,jj)=FidNoise(kk,jj)+abs((Operator*kron(Basis1(:,oo),Basis2(:,nn)))'*(OperId64*kron(Basis1(:,oo),Basis2(:,nn))))^2/16;
                        end
                    end
                    
                end
            end
            FidNoiseAv(ll)=sum(sum(FidNoise))/Edcnoisedim^2
        else
            ttgate(ll) = NaN;
            Fid16(ll) = NaN;
            FidNoiseAv(ll) = NaN;
        end
        save 2QFid_noise_adiabatic_Vt_05.mat
        toc;
        ll
%         (ll-1)/Vtdim+(mm-1)/Vtdim/EdcFdim+1/Vtdim/EdcFdim
end
profile viewer

%%
figure(67381);
subplot(2,2,1);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,ttgate(:,1,3)*1e9);
ylabel('Gate time (ns)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,2);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid16(:,1,3)))/log(10));
ylabel('log_{10}(Leakage error)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,3);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-FidNoiseAv(:,1,3))/log(10));
ylabel('log_{10}(Noise error)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,4);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(RelErr(:,1,3))/log(10));
ylabel('log_{10}(Relaxation error)')
xlabel('V_t - \epsilon_{ff} (MHz)')

%% Origin
VVt=Vt'/1e9;
Gtime=ttgate*1e9;
LeakN=1-Fid16;
TotN=1-FidNoiseAv;
RRelErr=RelErr;

%% sweep Vt, K, Edc and RMSnoise and calculate fidelity - Edc1=Edc2
profile on
clearvars

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C
e0=8.854187817e-12; %dielectric constant

% Energy unit is Hz (E/h)
A=117e6; %contact hyperfine
B0=0.2*2;
L1=7.5e-9; %half interdot distance
L2=7.5e-9; %half interdot distance
dg1=-0.002;%relative change in g-factor (delta_g/g)
dg2=-0.002;%relative change in g-factor (delta_g/g)
z=180e-9; %INTERDONOR SEPARATION
Edc0=-20000;

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0];
sigma_z=[1 0; 0 -1];

% define Hamiltonian matrices %%
% basis => orbital_1 x electron_1 x nucleus_1 x orbital_2 x electron_2 x nucleus_2  

H_Hyper1=(tensor(eye(2),sigma_x,sigma_x,eye(8))+tensor(eye(2),sigma_y,sigma_y,eye(8))+tensor(eye(2),sigma_z,sigma_z,eye(8)))*A/4*tensor(eye(2)/2-sigma_x/2,eye(32));
H_Hyper2=(tensor(eye(16),sigma_x,sigma_x)+tensor(eye(16),sigma_y,sigma_y)+tensor(eye(16),sigma_z,sigma_z))*A/4*tensor(eye(8),eye(2)/2-sigma_x/2,eye(4));
H_Znuc1=-tensor(eye(4),17.2e6*B0/2*sigma_z,eye(8));
H_Znuc2=-tensor(eye(32),17.2e6*B0/2*sigma_z);
H_Zeeman1=tensor(eye(2)+(eye(2)/2+sigma_x/2)*dg1,28e9*B0/2*sigma_z,eye(16));
H_Zeeman2=tensor(eye(8),eye(2)+(eye(2)/2+sigma_x/2)*dg2,28e9*B0/2*sigma_z,eye(2));
H_dipol=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2)*(tensor(sigma_x,eye(4),sigma_x,eye(4))+tensor(sigma_x,eye(32))+tensor(eye(8),sigma_x,eye(4)));%+tensor(identity,identity,identity,sigma_x,identity,identity)+tensor(identity,identity,sigma_x,identity,identity,identity)

Vtdim=10;
detmin=320e6;
detmax=460e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

EdcFdim=3;
EdcFmin=-400;
EdcFmax=0;
EdcF=(EdcFmin:(EdcFmax-EdcFmin)/(EdcFdim-1):EdcFmax);

Edcdim=2000; % Edc setup adiabatic sweep resolution

Edcnoisedim=3;

RMSdim=15;
RMSmin=10;
RMSmax=1000;
RMSnoise=logspace(log(RMSmin)/log(10),log(RMSmax)/log(10),RMSdim);

Kdim=15;
Kmin=1;
Kmax=100;
K=logspace(log(Kmin)/log(10),log(Kmax)/log(10),Kdim);

ttgate = zeros(EdcFdim,Vtdim,Kdim);
Fid16 = zeros(EdcFdim,Vtdim,Kdim,RMSdim,Edcnoisedim,Edcnoisedim);

dt=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);

gdd=1/(h*16*pi*e0*11.7*z^3)*(e*2*L1)*(e*2*L2);
dDorbdEz=e*2*L1*pi/h;

gff=zeros(1,Edcdim);

for ll = 1:Vtdim
    H_tunel1=tensor(Vt(ll)/2*sigma_z,eye(32));
    H_tunel2=tensor(eye(8),Vt(ll)/2*sigma_z,eye(4));
    H_0=H_Znuc1+H_Znuc2+H_tunel1+H_tunel2+H_Zeeman1+H_Zeeman2+H_Hyper1+H_Hyper2+H_dipol;
    gorb=Vt(ll)*pi;
    
    for mm = 1:EdcFdim
        tic;
        % find adiabatic t(E)
        dEdcMeas=(EdcF(mm)-Edc0)/(Edcdim-1);
        Edc=Edc0:dEdcMeas:EdcF(mm);
        
        for kk=1:Edcdim
            eo=sqrt(Vt(ll)^2+(e*Edc(kk)*2*L1/h+2*gdd)^2);
            pos=-(e*Edc(kk)*2*L1/h+2*gdd)/eo;
            eff=sqrt(((28e9*(1+dg1*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
            Dso(kk)=(eo-gdd*Vt(ll)/eo-eff)*pi;
            gso(kk)=A/4/sqrt(2)*Vt(ll)/eo*2*pi;
            Dorb=(e*Edc(kk)*2*L1/h+2*gdd)*pi; %IS THAT TRUE?
            dtOrb(kk)=1*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5; %IS THAT TRUE?
        end
        
        dDsodEz=-derivative(Dso)/dEdcMeas;
        dgsodEz=derivative(gso)/dEdcMeas;
        
        for kk=1:Edcdim
            dtSO(kk)=1*abs(dDsodEz(kk)*dEdcMeas*gso(kk)-dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
            dt1(kk)=max(dtOrb(kk),dtSO(kk));
        end
        
        % measure gff
        for kk=1:Edcdim
            H_Edc1=tensor(e*Edc(kk)*L1/h*sigma_x,eye(32));
            H_Edc2=tensor(eye(8),e*Edc(kk)*L2/h*sigma_x,eye(4));
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc(kk)*L1/h*sigma_x);
            singlet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])-tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
            triplet=1/sqrt(2)*(tensor(Orb_eVec(:,1)',[1 0],[0 1],Orb_eVec(:,1)',[0 1],[1 0])+tensor(Orb_eVec(:,1)',[0 1],[1 0],Orb_eVec(:,1)',[1 0],[0 1]));
            H=H_0+H_Edc1+H_Edc2;
            [H_eVec,H_eVal]=eig(H);
            proj1=0;
            proj2=0;
            aa=0;
            bb=0;
            for i=1:64
                AA=abs(singlet*H_eVec(:,i))^2;
                BB=abs(triplet*H_eVec(:,i))^2;
                if AA>proj1
                    proj1=AA;
                    aa=i;
                end
                if BB>proj2
                    proj2=BB;
                    bb=i;
                end
            end
            gff(kk)=abs((H_eVal(bb,bb)-H_eVal(aa,aa))/2);
        end
        
        for kkkk=1:Kdim
            dt=K(kkkk)*dt1;
            SetupRot=sum(gff.*dt)*4*pi;
            tgate=(pi/2-2*SetupRot)/4/pi/gff(Edcdim);
            
            if sign(tgate)>=0
                % 2-q gate time evolution
                ttgate(mm,ll,kkkk) = tgate + 2*sum(dt);
%                 ttgate(mm,ll,kkkk)
                
                [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc0*L1/h*sigma_x);
                UP1=tensor(Orb_eVec(:,1)',[1 0],[0 1])';
                DOWN1=tensor(Orb_eVec(:,1)',[0 1],[1 0])';
                UP2=UP1;
                DOWN2=DOWN1;
                
                OperStep=zeros(64,64,Edcdim);
                Operator=eye(64);
                
                for t=1:Edcdim
                    H_Edc1=tensor(e*Edc(t)*L1/h*sigma_x,eye(32));
                    H_Edc2=tensor(eye(8),e*Edc(t)*L2/h*sigma_x,eye(4));
                    H=H_0+H_Edc1+H_Edc2;
                    OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                    Operator=OperStep(:,:,t)*Operator;
                end
                
                Operator=expm(-1i*2*pi*H*tgate)*Operator;
                
                for t=(Edcdim+1):(2*Edcdim)
                    Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                end
                
                % fidelities
                Basis1(:,1)=UP1;
                Basis1(:,2)=DOWN1;
                Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
                Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
                Basis2(:,1)=UP2;
                Basis2(:,2)=DOWN2;
                Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
                Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
                Basis2q(:,1)=kron(UP1,UP2);
                Basis2q(:,2)=kron(UP1,DOWN2);
                Basis2q(:,3)=kron(DOWN1,UP2);
                Basis2q(:,4)=kron(DOWN1,DOWN2);
                
                phase=angle(Basis2q(:,2)'*Operator*Basis2q(:,2))-angle(Basis2q(:,4)'*Operator*Basis2q(:,4));
                
                OperId=expm(-1i*(-phase/2*(kron(eye(2),sigma_z)+kron(sigma_z,eye(2)))))*[1 0 0 0; 0 1/sqrt(2) -1i/sqrt(2) 0; 0 -1i/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];%expm(-1i*pi/4*tensor(sigma_x,sigma_x));
                OperId64=zeros(64);
                for ii=1:4
                    for jj=1:4
                        OperId64=OperId64+OperId(ii,jj)*Basis2q(:,ii)*Basis2q(:,jj)';
                    end
                end
                
                for ii=1:4
                    for jj=1:4
                        Fid16(mm,ll,kkkk,:,(Edcnoisedim+1)/2,(Edcnoisedim+1)/2)=Fid16(mm,ll,kkkk,:,(Edcnoisedim+1)/2,(Edcnoisedim+1)/2)+abs((Operator*kron(Basis1(:,ii),Basis2(:,jj)))'*(OperId64*kron(Basis1(:,ii),Basis2(:,jj))))^2/16;
                    end
                end
%                 Fid16(mm,ll,kkkk,:,(Edcnoisedim+1)/2,(Edcnoisedim+1)/2)
                
                % sweep Edc noise and average
                for rms=1:RMSdim
                    Edcnoise=-1:2/(Edcnoisedim-1):1;
                    Edcnoisemax=RMSnoise(rms)/sqrt(sum(Edcnoise.^2)/Edcnoisedim);
                    Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;
                    for nnn1=1:Edcnoisedim
                        Edc1=Edc+Edcnoise(nnn1);
                        [Orb_eVec1,Orb_eVal1]=eig(Vt(ll)/2*sigma_z+e*Edc1(1)*L1/h*sigma_x);
                        UP1=tensor(Orb_eVec1(:,1)',[1 0],[0 1])';
                        DOWN1=tensor(Orb_eVec1(:,1)',[0 1],[1 0])';
                        for nnn2=1:Edcnoisedim
                            if (nnn1~=(Edcnoisedim+1)/2)||(nnn2~=(Edcnoisedim+1)/2)
                                Edc2=Edc+Edcnoise(nnn2);
                                [Orb_eVec2,Orb_eVal2]=eig(Vt(ll)/2*sigma_z+e*Edc2(1)*L2/h*sigma_x);
                                UP2=tensor(Orb_eVec2(:,1)',[1 0],[0 1])';
                                DOWN2=tensor(Orb_eVec2(:,1)',[0 1],[1 0])';
                                
                                % time evolution
                                OperStep=zeros(64,64,Edcdim);
                                Operator=eye(64);
                                
                                for t=1:Edcdim
                                    H_Edc1=tensor(e*Edc1(t)*L1/h*sigma_x,eye(32));
                                    H_Edc2=tensor(eye(8),e*Edc2(t)*L2/h*sigma_x,eye(4));
                                    H=H_0+H_Edc1+H_Edc2;
                                    OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                                    Operator=OperStep(:,:,t)*Operator;
                                end
                                
                                Operator=expm(-1i*2*pi*H*tgate)*Operator;
                                
                                for t=(Edcdim+1):(2*Edcdim)
                                    Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
                                end
                                
                                % fidelities
                                Basis1(:,1)=UP1;
                                Basis1(:,2)=DOWN1;
                                Basis1(:,3)=(UP1+DOWN1)/sqrt(2);
                                Basis1(:,4)=(1i*UP1+DOWN1)/sqrt(2);
                                Basis2(:,1)=UP2;
                                Basis2(:,2)=DOWN2;
                                Basis2(:,3)=(UP2+DOWN2)/sqrt(2);
                                Basis2(:,4)=(1i*UP2+DOWN2)/sqrt(2);
                                Basis2q(:,1)=kron(UP1,UP2);
                                Basis2q(:,2)=kron(UP1,DOWN2);
                                Basis2q(:,3)=kron(DOWN1,UP2);
                                Basis2q(:,4)=kron(DOWN1,DOWN2);
                                
                                OperId64=zeros(64);
                                for oo=1:4
                                    for nn=1:4
                                        OperId64=OperId64+OperId(oo,nn)*Basis2q(:,oo)*Basis2q(:,nn)';
                                    end
                                end
                                
                                for oo=1:4
                                    for nn=1:4
                                        Fid16(mm,ll,kkkk,rms,nnn1,nnn2)=Fid16(mm,ll,kkkk,rms,nnn1,nnn2)+abs((Operator*kron(Basis1(:,oo),Basis2(:,nn)))'*(OperId64*kron(Basis1(:,oo),Basis2(:,nn))))^2/16;
                                    end
                                end
                                %                             Fid16(mm,ll,kkkk,rms,nnn1,nnn2)
                            end
                        end
                    end
                end
            else
                ttgate(mm,ll,kkkk) = NaN;
                Fid16(mm,ll,kkkk,:,:,:) = NaN;
            end
            save 2QFid_noise_adiabatic_RMS_K_09.mat
        end
%         (ll-1)/Vtdim+(mm-1)/Vtdim/EdcFdim+1/Vtdim/EdcFdim
        toc;
    end
end
profile viewer


%% plot gate time and error for an all K and RMSnoise values (only works if RMSdim=Kdim=3)
figure(6725311);
subplot(3,5,1);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,ttgate(:,:,1)'*1e9);axis xy;colorbar;
title('time (ns), K=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,2);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid16(:,:,1,1,(Edcnoisedim+1)/2,(Edcnoisedim+1)/2)))'/log(10));axis xy;colorbar;
title('log_{10}(leakage error), K=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,3);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,1,1,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=10, RMS=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,4);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,1,2,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=10, RMS=100')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,1,3,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=10, RMS=1000')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,1+5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,ttgate(:,:,2)'*1e9);axis xy;colorbar;
title('time (ns), K=45')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,2+5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid16(:,:,2,1,(Edcnoisedim+1)/2,(Edcnoisedim+1)/2)))'/log(10));axis xy;colorbar;
title('log_{10}(leakage error), K=45')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,3+5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,2,1,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=45, RMS=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,4+5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,2,2,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=45, RMS=100')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,5+5);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,2,3,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=45, RMS=1000')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,1+10);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,ttgate(:,:,3)'*1e9);axis xy;colorbar;
title('time (ns), K=200')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,2+10);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid16(:,:,3,1,(Edcnoisedim+1)/2,(Edcnoisedim+1)/2)))'/log(10));axis xy;colorbar;
title('log_{10}(leakage error), K=200')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,3+10);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,3,1,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=200, RMS=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,4+10);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,3,2,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=200, RMS=100')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,5+10);
imagesc(-EdcF,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(sum(Fid16(:,:,3,3,:,:),6),5)/Edcnoisedim^2)'/log(10));axis xy;colorbar;
title('log_{10}(total error), K=200, RMS=1000')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
 
%% plot optimized total error vs RMSnoise and K - optimze Edc as well
%figure 6c
%tpi2 = zeros(Edcdim,Vtdim,Kdim);
%Fid4 = zeros(Edcdim,Vtdim,Kdim,RMSdim,Edcnoisedim);
TotErr=zeros(Kdim,RMSdim);
VVt=zeros(Kdim,RMSdim);
EEdc=zeros(Kdim,RMSdim);
Gtime=zeros(Kdim,RMSdim);
 
for kkk=1:Kdim
    for rms=1:RMSdim
        [mm,eedc]=min(1-sum(sum(Fid16(:,:,kkk,rms,:,:),6),5)/Edcnoisedim^2);
        [mm,vvt]=min(mm);
        eedc=eedc(vvt);
        EEdc(kkk,rms)=-EdcF(eedc);
        VVt(kkk,rms)=Vt(vvt);
        Gtime(kkk,rms)=ttgate(eedc,vvt,kkk);
        TotErr(kkk,rms)=mm;
    end
end
    
figure(67253211);
subplot(2,2,1);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),log(TotErr)/log(10));axis xy;colorbar;
title('total error')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')
subplot(2,2,2);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),Gtime);axis xy;colorbar;
title('optimal gate time (ns)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')
subplot(2,2,3);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),(VVt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6);axis xy;colorbar;
title('optimal Vt-eff (MHz)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')
subplot(2,2,4);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),EEdc);axis xy;colorbar;
title('optimal EdcF (V/m)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')

%% Origin
Rnoise=RMSnoise;
Kk=K';
TTotErr=TotErr;
GGtime=sum(Gtime,2)/RMSdim;
min(min(TTotErr))
max(max(TTotErr))
 
%% plot optimized total error vs RMSnoise and K - assume Edc=0 - MANUALLY SELECT EDC=0 INDEX!!!!!!!!!!!!!
%tpi2 = zeros(Edcdim,Vtdim,Kdim);
%Fid4 = zeros(Edcdim,Vtdim,Kdim,RMSdim,Edcnoisedim);
TotErr=zeros(Kdim,RMSdim);
VVt=zeros(Kdim,RMSdim);
EEdc=zeros(Kdim,RMSdim);
Gtime=zeros(Kdim,RMSdim);
 
eedc=2; % MANUALLY SELECT EDC=0 INDEX!!!!!!!!!!!!!
 
for kkk=1:Kdim
    for rms=1:RMSdim
        [mm,vvt]=min(1-sum(sum(Fid16(eedc,:,kkk,rms,:,:),6),5)/Edcnoisedim^2);
        EEdc(kkk,rms)=-EdcF(eedc);
        VVt(kkk,rms)=Vt(vvt);
        Gtime(kkk,rms)=ttgate(eedc,vvt,kkk);
        TotErr(kkk,rms)=mm;
    end
end
    
figure(672532111);
subplot(2,2,1);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),log(TotErr)/log(10));axis xy;colorbar;
title('total error')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')
subplot(2,2,2);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),Gtime);axis xy;colorbar;
title('optimal gate time (ns)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')
subplot(2,2,3);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),VVt);axis xy;colorbar;
title('optimal Vt (GHz)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')
subplot(2,2,4);
imagesc(log(RMSnoise)/log(10),log(K)/log(10),EEdc);axis xy;colorbar;
title('optimal EdcF (V/m)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')



