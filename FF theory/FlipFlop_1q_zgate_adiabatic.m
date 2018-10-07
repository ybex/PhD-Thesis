%% setup parameters
clear all;
addpath Q:\spin-QED\Theory_matlab\Matlab_functions
addpath Q:\spin-QED\Papers\Quantum processor\NatComms\Matlab simulations

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C 

A=117e6;
B0=0.2*2;
Vt=11.4331760406494e9;%sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+173e6; %double-dot tunel rate
Edc0=-125000*B0;%-110;%-330;
L=7.5e-9; %half interdot distance
dg=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% basis =>  nucleus x orbital x electron
% 01) u,e,u
% 02) u,e,d
% 03) u,g,u
% 04) u,g,d
% 05) d,e,u
% 06) d,e,d
% 07) d,g,u
% 08) d,g,d

% define Hamiltonian matrices %%
H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);
H_0=H_Znuc+H_tunel+H_Zel+H_Hyper;

%% find adiabatic t(E)
EdcF=-403;
Edcdim=50000;
dEdcMeas=(EdcF-Edc0)/(Edcdim-1);
Edc=Edc0:dEdcMeas:EdcF;

dt=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);
K=50;

dDorbdEz=e*2*L*pi/h;
gorb=Vt*pi;

for kk=1:Edcdim
    eo=sqrt(Vt^2+(e*Edc(kk)*2*L/h)^2);
    pos=-(e*Edc(kk)*2*L/h)/eo;
    eff=sqrt(((28e9*(1+dg*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
    Dso(kk)=(eo-eff)*pi;
    gso(kk)=A/4*Vt/eo*2*pi;
    Dorb=(e*Edc(kk)*2*L/h)*pi;
    dtOrb(kk)=K*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5;
end

dDsodEz=-derivative(Dso)/dEdcMeas;
dgsodEz=derivative(gso)/dEdcMeas;

for kk=1:Edcdim
    dtSO(kk)=K*(dDsodEz(kk)*dEdcMeas*gso(kk)-dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
    dt(kk)=max(dtOrb(kk),dtSO(kk));
end

figure(1123);
subplot(1,2,1);
semilogy(Edc,dt./dEdcMeas,Edc,dtSO./dEdcMeas,Edc,dtOrb./dEdcMeas)
legend('dt','dtSO','dtOrb')
ylabel('dt/dEz (s)')
xlabel('Edc (V/m)')

%% time-evolution

gatedim=5000;

%initial basis
H_Edc=tensor(identity,e*Edc0*L/h*sigma_x,identity);
H=H_0+H_Edc;
[H_eVec,H_eVal]=eig(H);
[Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc0*L/h*sigma_x);
UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);
DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);
proj1=0;
proj2=0;
aa=0;
bb=0;
for i=1:8
    AA=abs(UP*H_eVec(:,i))^2;
    BB=abs(DOWN*H_eVec(:,i))^2;
    if AA>proj1
        proj1=AA;
        aa=i;
    end
    if BB>proj2
        proj2=BB;
        bb=i;
    end
end
UP=sign(UP*H_eVec(:,aa))*H_eVec(:,aa);
DOWN=sign(DOWN*H_eVec(:,bb))*H_eVec(:,bb);

OperStep=zeros(8,8,Edcdim);
Operator=eye(8);
nuE=zeros(1,Edcdim);

down=tensor([1 0]'*[1 0],eye(2),[0 1]'*[0 1]);
up=tensor([0 1]'*[0 1],eye(2),[1 0]'*[1 0]);
int=tensor(eye(2),[1 -1]'*[1 -1]./2,eye(2));
don=tensor(eye(2),[1 1]'*[1 1]./2,eye(2));
charup=tensor(eye(2),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(2));
%chardown=tensor(eye(2),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(2));

Psit = zeros(8,2*Edcdim+gatedim+1);
Psit(:,1)=(1*UP+1*DOWN)/sqrt(2);

FFqubit=zeros(1,2*Edcdim+gatedim+1);
FFx=zeros(1,2*Edcdim+gatedim+1);
position=zeros(1,2*Edcdim+gatedim+1);
ge=zeros(1,2*Edcdim+gatedim+1);
ge(1)=Psit(:,1)'*(charup)*Psit(:,1);
FFqubit(1)=Psit(:,1)'*(up-down)*Psit(:,1);
position(1)=Psit(:,1)'*(don-int)*Psit(:,1);
xplus=(tensor([0 1],[1 0],[1 0])+tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+tensor([1 0],[0 1],[0 1]))/sqrt(2);
xminus=(tensor([0 1],[1 0],[1 0])-tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])-tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])-tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])-tensor([1 0],[0 1],[0 1]))/sqrt(2);
FFx(1)=Psit(:,1)'*(xplus-xminus)*Psit(:,1);

ttt=zeros(1,2*Edcdim+gatedim+1);
for t=1:Edcdim
    H_Edc=tensor(identity,e*Edc(t)*L/h*sigma_x,identity);
    H=H_0+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc(t)*L/h*sigma_x);
    charup=tensor(eye(2),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(2));
    %chardown=tensor(eye(2),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(2));
    UP1=tensor([0 1],Orb_eVec(:,1)',[1 0]);%H_eVec(:,3);%
    DOWN1=tensor([1 0],Orb_eVec(:,1)',[0 1]);%H_eVec(:,1);%
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:8
        AA=abs(UP1*H_eVec(:,i))^2;
        BB=abs(DOWN1*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    nuE(t)=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
    OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
    Operator=OperStep(:,:,t)*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*(charup)*Psit(:,1+t);
    FFqubit(1+t)=Psit(:,1+t)'*(up-down)*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*(don-int)*Psit(:,1+t);
    ttt(t+1)=ttt(t)+dt(t);
    xplus=(tensor([0 1],[1 0],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[0 1],[0 1]))/sqrt(2);
    xminus=(tensor([0 1],[1 0],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[0 1],[0 1]))/sqrt(2);
    FFx(1+t)=Psit(:,1+t)'*(xplus-xminus)*Psit(:,1+t);
end

SetupRot=sum((nuE-nuE(1)).*dt)*2*pi
%figure(5648)
%subplot(2,1,1)
%plot(Edc,nuE-nuE(1))

Gate=1*pi;
tgate=(Gate/2/pi-SetupRot/pi)/(nuE(Edcdim)-nuE(1))%pi z-gate
OperStepp=expm(-1i*2*pi*H*tgate/gatedim);

for t=(Edcdim+1):(Edcdim+gatedim)
    Operator=OperStepp*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*(charup)*Psit(:,1+t);
    FFqubit(1+t)=Psit(:,1+t)'*(up-down)*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*(don-int)*Psit(:,1+t);
    ttt(t+1)=ttt(t)+tgate/gatedim;
    xplus=(tensor([0 1],[1 0],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[0 1],[0 1]))/sqrt(2);
    xminus=(tensor([0 1],[1 0],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[0 1],[0 1]))/sqrt(2);
    FFx(1+t)=Psit(:,1+t)'*(xplus-xminus)*Psit(:,1+t);
end

for t=(Edcdim+gatedim+1):(2*Edcdim+gatedim)
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc(2*Edcdim+gatedim-t+1)*L/h*sigma_x);
    charup=tensor(eye(2),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(2));
    %chardown=tensor(eye(2),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(2));
    Operator=OperStep(:,:,2*Edcdim+gatedim-t+1)*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*(charup)*Psit(:,1+t);
    FFqubit(1+t)=Psit(:,1+t)'*(up-down)*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*(don-int)*Psit(:,1+t);
    ttt(t+1)=ttt(t)+dt(2*Edcdim+gatedim-t+1);
    xplus=(tensor([0 1],[1 0],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+exp(-1i*2*pi*nuE(1)*ttt(t+1))*tensor([1 0],[0 1],[0 1]))/sqrt(2);
    xminus=(tensor([0 1],[1 0],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[1 0],[0 1]))'/sqrt(2)*(tensor([0 1],[1 0],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[1 0],[0 1]))/sqrt(2)+(tensor([0 1],[0 1],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[0 1],[0 1]))'/sqrt(2)*(tensor([0 1],[0 1],[1 0])+exp(-1i*(2*pi*nuE(1)*ttt(t+1)+pi))*tensor([1 0],[0 1],[0 1]))/sqrt(2);
    FFx(1+t)=Psit(:,1+t)'*(xplus-xminus)*Psit(:,1+t);
end

%%
Edcgate=zeros(1,gatedim);
Edcgate(:)=Edc(Edcdim);

FigHandle = figure(5648);
set(FigHandle, 'Position', [500, 100, 900, 900]);

subplot(5,1,1);
plot(ttt*1e9,[Edc Edcgate Edc(Edcdim:-1:1) Edc(1)]/1e3,'LineWidth',2);%,'Color',[.5 .5 .5]
ylabel('Ez-Ez0 (kV/m)')
xlim([-1 69.7])
ylim([-55 5])
set(gca,'YTick',[-40 0])%,'Position', [0.12 0.75 0.8 0.2]

subplot(5,1,2);
plot(ttt*1e9,FFqubit);
ylabel('FFz')
xlim([-1 69.7])
ylim([-0.015 0.009])
set(gca,'YTick',[-0.01 0])

subplot(5,1,3);
plot(ttt*1e9,FFx,'LineWidth',2);
ylabel('FFx_{rot}')
xlim([-1 69.7])
ylim([-1.2 1.2])
set(gca,'YTick',[-1 1])

subplot(5,1,4);
plot(ttt*1e9,position);
ylabel('position')
xlim([-1 69.7])
ylim([-0.1 1.1])
set(gca,'YTick',[0 1])

subplot(5,1,5);
plot(ttt*1e9,ge);
xlabel('Time (ns)')
ylabel('Charge')
xlim([-1 69.7])
ylim([-0.0005 0.0055])
set(gca,'YTick',[0 0.005])

%% fidelitiy

Basis(:,1)=UP;
Basis(:,2)=DOWN;
Basis(:,3)=(UP+DOWN)/sqrt(2);
Basis(:,4)=(1i*UP+DOWN)/sqrt(2);

phase=Gate+nuE(1)*ttt(2*Edcdim+gatedim+1)*2*pi;
OperId=expm(-1i*phase/2*sigma_z);
OperId8=zeros(8,8);
for ii=1:2
    for jj=1:2
        OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
    end
end

Fid4 = 0;
for ii=1:4
        Fid4=Fid4+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
end
Error=1-Fid4

%% sweep Edc noise and average        
Edcnoisedim=3;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=100/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;
Error4 = zeros(1,Edcnoisedim);

for kk=1:Edcnoisedim
    kk
    H_Edc=tensor(identity,e*(Edc0+Edcnoise(kk))*L/h*sigma_x,identity);
    H=H_0+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*(Edc0+Edcnoise(kk))*L/h*sigma_x);
    UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);
    DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:8
        AA=abs(UP*H_eVec(:,i))^2;
        BB=abs(DOWN*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    UP=sign(UP*H_eVec(:,aa))*H_eVec(:,aa);
    DOWN=sign(DOWN*H_eVec(:,bb))*H_eVec(:,bb);

    OperStep=zeros(8,8,Edcdim);
    Operator=eye(8);
    
    for t=1:Edcdim
        H_Edc=tensor(identity,e*(Edc(t)+Edcnoise(kk))*L/h*sigma_x,identity);
        H=H_0+H_Edc;
        OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
        Operator=OperStep(:,:,t)*Operator;
    end
    
    OperStepp=expm(-1i*2*pi*H*tgate);
    Operator=OperStepp*Operator;
    
    for t=(Edcdim+gatedim+1):(2*Edcdim+gatedim)
        Operator=OperStep(:,:,2*Edcdim+gatedim-t+1)*Operator;
    end

    % fidelities
    Basis(:,1)=UP;
    Basis(:,2)=DOWN;
    Basis(:,3)=(UP+DOWN)/sqrt(2);
    Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
    
    OperId8=zeros(8,8);
    for ii=1:2
        for jj=1:2
            OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
        end
    end
    
    Error4(kk) = 0;
    for ii=1:4
        Error4(kk)=Error4(kk)+1/4-abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
    end
    Error4(kk)
end
ErrorAv=sum(Error4)/Edcnoisedim

%% calculate fidelity and gate time as a function of desired phase gate
Gatedim=20;
Gate=0:2*pi/(Gatedim-1):2*pi;

Edcnoisedim=15;
Edcnoisemax=170;
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

timegate=zeros(1,Gatedim);
phasegate=zeros(1,Gatedim);
Error=zeros(Gatedim,Edcnoisedim);

OperStep=zeros(8,8,Edcdim,Edcnoisedim);
OperatorSetup=zeros(8,8,Edcnoisedim);
Operator=zeros(8,8,Edcnoisedim);

%initial basis
H_Edc=tensor(identity,e*Edc0*L/h*sigma_x,identity);
H=H_0+H_Edc;
[H_eVec,H_eVal]=eig(H);
[Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc0*L/h*sigma_x);
UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);
DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);
proj1=0;
proj2=0;
aa=0;
bb=0;
for i=1:8
    AA=abs(UP*H_eVec(:,i))^2;
    BB=abs(DOWN*H_eVec(:,i))^2;
    if AA>proj1
        proj1=AA;
        aa=i;
    end
    if BB>proj2
        proj2=BB;
        bb=i;
    end
end
UP=sign(UP*H_eVec(:,aa))*H_eVec(:,aa);
DOWN=sign(DOWN*H_eVec(:,bb))*H_eVec(:,bb);
Basis(:,1)=UP;
Basis(:,2)=DOWN;
Basis(:,3)=(UP+DOWN)/sqrt(2);
Basis(:,4)=(1i*UP+DOWN)/sqrt(2);

kk=(Edcnoisedim+1)/2

OperatorSetup(:,:,kk)=eye(8);

for t=1:Edcdim
    H_Edc=tensor(identity,e*Edc(t)*L/h*sigma_x,identity);
    H=H_0+H_Edc;
    OperStep(:,:,t,kk)=expm(-1i*2*pi*H*dt(t));
    OperatorSetup(:,:,kk)=OperStep(:,:,t,kk)*OperatorSetup(:,:,kk);
end

for gg=1:Gatedim
    tgate=Gate(gg)/2/pi/(nuE(Edcdim)-nuE(1));
    Operator(:,:,kk)=expm(-1i*2*pi*H*tgate)*OperatorSetup(:,:,kk);
    for t=(Edcdim+1):(2*Edcdim)
        Operator(:,:,kk)=OperStep(:,:,2*Edcdim-t+1,kk)*Operator(:,:,kk);
    end
    
    % fidelitiy
    setuptime=0;
    for eee=1:Edcdim
        setuptime=setuptime+dt(eee);
    end
    %SetupRot=sum((nuE-nuE(1)).*dt)*2*pi
    phasegate(gg)=Gate(gg)+2*SetupRot;
    timegate(gg)=tgate+2*setuptime;
    phase=phasegate(gg)+nuE(1)*(timegate(gg))*2*pi;
    
    OperId=expm(-1i*phase/2*sigma_z);
    OperId8=zeros(8,8);
    for ii=1:2
        for jj=1:2
            OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
        end
    end
    
    Fid4 = 0;
    for ii=1:4
        Fid4=Fid4+abs((Operator(:,:,kk)*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
    end
    Error(gg,kk)=1-Fid4;
end

% sweep Edc noise and average    
for kk=[1:(Edcnoisedim-1)/2 (Edcnoisedim+3)/2:Edcnoisedim]
    kk
    H_Edc=tensor(identity,e*(Edc0+Edcnoise(kk))*L/h*sigma_x,identity);
    H=H_0+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*(Edc0+Edcnoise(kk))*L/h*sigma_x);
    UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);
    DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:8
        AA=abs(UP*H_eVec(:,i))^2;
        BB=abs(DOWN*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    UP=sign(UP*H_eVec(:,aa))*H_eVec(:,aa);
    DOWN=sign(DOWN*H_eVec(:,bb))*H_eVec(:,bb);
    Basis(:,1)=UP;
    Basis(:,2)=DOWN;
    Basis(:,3)=(UP+DOWN)/sqrt(2);
    Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
    
    OperatorSetup(:,:,kk)=eye(8);
    
    for t=1:Edcdim
        H_Edc=tensor(identity,e*(Edc(t)+Edcnoise(kk))*L/h*sigma_x,identity);
        H=H_0+H_Edc;
        OperStep(:,:,t,kk)=expm(-1i*2*pi*H*dt(t));
        OperatorSetup(:,:,kk)=OperStep(:,:,t,kk)*OperatorSetup(:,:,kk);
    end
    
    for gg=1:Gatedim
        tgate=Gate(gg)/2/pi/(nuE(Edcdim)-nuE(1));
        Operator(:,:,kk)=expm(-1i*2*pi*H*tgate)*OperatorSetup(:,:,kk);
        for t=(Edcdim+1):(2*Edcdim)
            Operator(:,:,kk)=OperStep(:,:,2*Edcdim-t+1,kk)*Operator(:,:,kk);
        end
        
        % fidelitiy
        phase=phasegate(gg)+nuE(1)*(timegate(gg))*2*pi;
        OperId=expm(-1i*phase/2*sigma_z);
        OperId8=zeros(8,8);
        for ii=1:2
            for jj=1:2
                OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
            end
        end
        
        Fid4 = 0;
        for ii=1:4
            Fid4=Fid4+abs((Operator(:,:,kk)*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
        end
        Error(gg,kk)=1-Fid4;
    end
end
ErrorAv=sum(Error,2)/Edcnoisedim;

% timegate
% phasegate
Error(:,(Edcnoisedim-1)/2).'
% ErrorAv.'
figure(0987)
plotyy(phasegate,timegate*1e9,phasegate,ErrorAv)


%% Calculate averaged error of a Setup+pi gate as a function of adiabaticity K and Edcnoisemax
Kdim=31;
Kmin=1;
Kmax=500;
K=logspace(log(Kmin)/log(10),log(Kmax)/log(10),Kdim);

Gate=pi;

Noisedim=21;
Noisemin=10;
Noisemax=1000;
Noise=logspace(log(Noisemin)/log(10),log(Noisemax)/log(10),Noisedim);
Edcnoisedim=3;

timegate=zeros(1,Kdim);
Error=zeros(Kdim,Noisedim,Edcnoisedim);

OperStep=zeros(8,8,Edcdim);

Edcdim=100000;
dEdcMeas=(EdcF-Edc0)/(Edcdim-1);
Edc=Edc0:dEdcMeas:EdcF;

% find adiabatic t(E) for K=1
dt1=zeros(1,Edcdim);
dtOrb=zeros(1,Edcdim);
dtSO=zeros(1,Edcdim);
Dso=zeros(1,Edcdim);
gso=zeros(1,Edcdim);

dDorbdEz=e*2*L*pi/h;
gorb=Vt*pi;

for kk=1:Edcdim
    eo=sqrt(Vt^2+(e*Edc(kk)*2*L/h)^2);
    pos=-(e*Edc(kk)*2*L/h)/eo;
    eff=sqrt(((28e9*(1+dg*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
    Dso(kk)=(eo-eff)*pi;
    gso(kk)=A/4*Vt/eo*2*pi;
    Dorb=(e*Edc(kk)*2*L/h)*pi;
	dtOrb(kk)=1*dDorbdEz*dEdcMeas*gorb/(Dorb^2+gorb^2)^1.5;
end
dDsodEz=-derivative(Dso)/dEdcMeas;
dgsodEz=derivative(gso)/dEdcMeas;
for kk=1:Edcdim
    dtSO(kk)=1*(dDsodEz(kk)*dEdcMeas*gso(kk)-dgsodEz(kk)*dEdcMeas*Dso(kk))/(Dso(kk)^2+gso(kk)^2)^1.5;
    dt1(kk)=max(dtOrb(kk),dtSO(kk));
end

setuptime1=0;
for eee=1:Edcdim
    setuptime1=setuptime1+dt1(eee);
end

% first without noise
tic;
kk=(Edcnoisedim+1)/2

nuE=zeros(1,Edcdim);

%initial basis
H_Edc=tensor(identity,e*Edc0*L/h*sigma_x,identity);
H=H_0+H_Edc;
[H_eVec,H_eVal]=eig(H);
[Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc0*L/h*sigma_x);
UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);
DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);
proj1=0;
proj2=0;
aa=0;
bb=0;
for i=1:8
    AA=abs(UP*H_eVec(:,i))^2;
    BB=abs(DOWN*H_eVec(:,i))^2;
    if AA>proj1
        proj1=AA;
        aa=i;
    end
    if BB>proj2
        proj2=BB;
        bb=i;
    end
end
nuE(1)=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
UP=sign(UP*H_eVec(:,aa))*H_eVec(:,aa);
DOWN=sign(DOWN*H_eVec(:,bb))*H_eVec(:,bb);
Basis(:,1)=UP;
Basis(:,2)=DOWN;
Basis(:,3)=(UP+DOWN)/sqrt(2);
Basis(:,4)=(1i*UP+DOWN)/sqrt(2);

%find nuE(Edc)
for t=2:Edcdim
    H_Edc=tensor(identity,e*Edc(t)*L/h*sigma_x,identity);
    H=H_0+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc(t)*L/h*sigma_x);
    UP1=tensor([0 1],Orb_eVec(:,1)',[1 0]);%H_eVec(:,3);%
    DOWN1=tensor([1 0],Orb_eVec(:,1)',[0 1]);%H_eVec(:,1);%
    proj1=0;
    proj2=0;
    aa=0;
    bb=0;
    for i=1:8
        AA=abs(UP1*H_eVec(:,i))^2;
        BB=abs(DOWN1*H_eVec(:,i))^2;
        if AA>proj1
            proj1=AA;
            aa=i;
        end
        if BB>proj2
            proj2=BB;
            bb=i;
        end
    end
    nuE(t)=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
end

tgate=Gate/2/pi/(nuE(Edcdim)-nuE(1));

for KK=1:Kdim
    setuptime=K(KK)*setuptime1;
    timegate(KK)=tgate+2*setuptime;
    dt=K(KK)*dt1;
    SetupRot=sum((nuE-nuE(1)).*dt)*2*pi;
    Operator=eye(8);
    for t=1:Edcdim
        H_Edc=tensor(identity,e*Edc(t)*L/h*sigma_x,identity);
        H=H_0+H_Edc;
        OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
        Operator=OperStep(:,:,t)*Operator;
    end
    Operator=expm(-1i*2*pi*H*tgate)*Operator;
    for t=(Edcdim+1):(2*Edcdim)
        Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
    end
    
    % fidelities
    phase=Gate+2*SetupRot+nuE(1)*timegate(KK)*2*pi;
    OperId=expm(-1i*phase/2*sigma_z);
    OperId8=zeros(8,8);
    for ii=1:2
        for jj=1:2
            OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
        end
    end
    
    Fid4 = 0;
    for ii=1:4
        Fid4=Fid4+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
    end
    Error(KK,:,kk)=1-Fid4;
end
%save 1Q-zgate_noise_adiabatic_09.mat
toc;
% sweep noise

for nn=1:Noisedim
    nn
    Edcnoise=-1:2/(Edcnoisedim-1):1;
    Edcnoisemax=Noise(nn)/sqrt(sum(Edcnoise.^2)/Edcnoisedim);
    Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;
    for kk=[1:(Edcnoisedim-1)/2 (Edcnoisedim+3)/2:Edcnoisedim]
        tic;
        kk
        %initial basis
        H_Edc=tensor(identity,e*(Edc0+Edcnoise(kk))*L/h*sigma_x,identity);
        H=H_0+H_Edc;
        [H_eVec,H_eVal]=eig(H);
        [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*(Edc0+Edcnoise(kk))*L/h*sigma_x);
        UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);
        DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);
        proj1=0;
        proj2=0;
        aa=0;
        bb=0;
        for i=1:8
            AA=abs(UP*H_eVec(:,i))^2;
            BB=abs(DOWN*H_eVec(:,i))^2;
            if AA>proj1
                proj1=AA;
                aa=i;
            end
            if BB>proj2
                proj2=BB;
                bb=i;
            end
        end
        UP=sign(UP*H_eVec(:,aa))*H_eVec(:,aa);
        DOWN=sign(DOWN*H_eVec(:,bb))*H_eVec(:,bb);
        Basis(:,1)=UP;
        Basis(:,2)=DOWN;
        Basis(:,3)=(UP+DOWN)/sqrt(2);
        Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
        
        for KK=1:Kdim
            dt=K(KK)*dt1;
            SetupRot=sum((nuE-nuE(1)).*dt)*2*pi;
            Operator=eye(8);
            for t=1:Edcdim
                H_Edc=tensor(identity,e*(Edc(t)+Edcnoise(kk))*L/h*sigma_x,identity);
                H=H_0+H_Edc;
                OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
                Operator=OperStep(:,:,t)*Operator;
            end
            Operator=expm(-1i*2*pi*H*tgate)*Operator;
            for t=(Edcdim+1):(2*Edcdim)
                Operator=OperStep(:,:,2*Edcdim-t+1)*Operator;
            end
            
            % fidelities
            phase=Gate+2*SetupRot+nuE(1)*timegate(KK)*2*pi;
            OperId=expm(-1i*phase/2*sigma_z);
            OperId8=zeros(8,8);
            for ii=1:2
                for jj=1:2
                    OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
                end
            end
            
            Fid4 = 0;
            for ii=1:4
                Fid4=Fid4+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
            end
            Error(KK,nn,kk)=1-Fid4;
        end
        save 1Q-zgate_noise_adiabatic_12.mat
        toc;
    end
end
%%
figure(1235)
subplot(1,2,1)
loglog(K,Error(:,1,(Edcnoisedim+1)/2))
xlabel('Average gate time (ns)')
ylabel('Error rate due to leakage')
ax1=gca;
XTickLabel=strings(1,Kdim);
for xtl=1:Kdim
XTickLabel(xtl)=num2str(timegate(xtl)*1e9);
end
set(ax1,'XTick',K,'XTickLabel',XTickLabel)
ax2 = axes('Position', get(ax1, 'Position'));
loglog(K,Error(:,1,(Edcnoisedim+1)/2))
set(ax2, 'XAxisLocation', 'top');
for xtl=1:Kdim
XTickLabel(xtl)=num2str(K(xtl));
end
set(ax2,'XTick',K,'XTickLabel',XTickLabel);
xlabel('adiabaticity K')

%RMSfactor=sqrt(sum((-1:2/(Edcnoisedim-1):1).^2)/Edcnoisedim);
subplot(1,2,2)
imagesc(log(Noise)/log(10),log(K)/log(10),log(sum(Error,3)/Edcnoisedim)/log(10)); axis xy; colorbar
%log(Noise*RMSfactor)
ylabel('log_{10}(K)')
xlabel('log_{10}(E_{noise}^{rms}) (V/m)')
title('log_{10}(Averaged error)')

%% origin
ErrorLeakage=Error(:,1,(Edcnoisedim+1)/2);
KKK=K';
SSetuptime=K'*setuptime1*1e9;
setuptime1*1e9

ErrorAvTable=sum(Error,3)/Edcnoisedim;
NNoise=Noise;%*RMSfactor
min(min(ErrorAvTable))
max(max(ErrorAvTable))

%% Charge qubit only - find adiabatic t(E)
clear all;
%addpath Z:\spin-QED\Theory_matlab\Matlab_functions

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C 

Vt=sqrt(((28e9+17.2e6-28e6)*1)^2 + 117^2/4)+173e6;%A*10/4;%-3e8;%6e9; %double-dot tunel rate
Edc0=-20000;%-110;%-330;
L=7.5e-9; %half interdot distance

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% define Hamiltonian matrices %%
H_tunnel=Vt/2*sigma_z;

% characterize adiabaticity
EdcF=20000;%-520;
Edcdim=10000;
dEdcMeas=(EdcF-Edc0)/(Edcdim-1);
EdcMeas=Edc0:dEdcMeas:EdcF;

dt=zeros(1,Edcdim);
K=100;

dDorbdEz=e*2*L*pi/h;
gorb=Vt*pi;
for kk=1:Edcdim
    Dorb=(e*EdcMeas(kk)*2*L/h)*pi;
    dt(kk)=K*dDorbdEz*gorb/(Dorb^2+gorb^2)^1.5;
end

figure(1121);
subplot(1,2,1);
plot(EdcMeas,dt)
ylabel('dt (s)')
xlabel('Edc (V/m)')

% time evolution
Edc=EdcMeas;
gatedim=1000;

%initial basis
H_Edc=e*Edc0*L/h*sigma_x;
H=H_tunnel+H_Edc;
[H_eVec,H_eVal]=eig(H);
UP=H_eVec(:,2)';
DOWN=H_eVec(:,1)';

OperStep=zeros(2,2,Edcdim);
Operator=eye(2);
nuE=zeros(1,Edcdim);

int=[1 -1]'*[1 -1]./2;
don=[1 1]'*[1 1]./2;
charup=H_eVec(:,2)*H_eVec(:,2)';
chardown=H_eVec(:,1)*H_eVec(:,1)';

Psit = zeros(2,2*Edcdim+gatedim+1);
Psit(:,1)=(0*UP+1*DOWN)/sqrt(1);

position=zeros(1,2*Edcdim+gatedim+1);
ge=zeros(1,2*Edcdim+gatedim+1);
ge(1)=Psit(:,1)'*charup*Psit(:,1)-Psit(:,1)'*chardown*Psit(:,1);
position(1)=Psit(:,1)'*don*Psit(:,1)-Psit(:,1)'*int*Psit(:,1);

for t=1:Edcdim
    H_Edc=e*Edc(t)*L/h*sigma_x;
    H=H_tunnel+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    charup=H_eVec(:,2)*H_eVec(:,2)';
    chardown=H_eVec(:,1)*H_eVec(:,1)';
    UP1=H_eVec(:,1)';
    DOWN1=H_eVec(:,1)';
    nuE(t)=abs(H_eVal(2,2)-H_eVal(1,1));
    OperStep(:,:,t)=expm(-1i*2*pi*H*dt(t));
    Operator=OperStep(:,:,t)*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*charup*Psit(:,1+t)-Psit(:,1+t)'*chardown*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*don*Psit(:,1+t)-Psit(:,1+t)'*int*Psit(:,1+t);
end

tgate=0.1e-9;
OperStepp=expm(-1i*2*pi*H*tgate/gatedim);

for t=(Edcdim+1):(Edcdim+gatedim)
    Operator=OperStepp*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*charup*Psit(:,1+t)-Psit(:,1+t)'*chardown*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*don*Psit(:,1+t)-Psit(:,1+t)'*int*Psit(:,1+t);
end

for t=(Edcdim+gatedim+1):(2*Edcdim+gatedim)
    [H_eVec,H_eVal]=eig(Vt/2*sigma_z+e*Edc(2*Edcdim+gatedim-t+1)*L/h*sigma_x);
    charup=H_eVec(:,2)*H_eVec(:,2)';
    chardown=H_eVec(:,1)*H_eVec(:,1)';
    Operator=OperStep(:,:,2*Edcdim+gatedim-t+1)*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*charup*Psit(:,1+t)-Psit(:,1+t)'*chardown*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*don*Psit(:,1+t)-Psit(:,1+t)'*int*Psit(:,1+t);
end

ttt=zeros(1,2*Edcdim+gatedim+1);
for kk=2:(Edcdim+1)
    ttt(kk)=ttt(kk-1)+dt(kk-1);
end

for kk=(Edcdim+2):(Edcdim+gatedim+1)
    ttt(kk)=ttt(kk-1)+tgate/gatedim;
end

for kk=(Edcdim+gatedim+2):(2*Edcdim+gatedim+1)
    ttt(kk)=ttt(kk-1)+dt(2*Edcdim+gatedim-kk+2);
end

figure(1121);
subplot(1,2,2)
plot(ttt,position,ttt,ge);
legend('position','Charge')