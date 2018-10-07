%% setup parameters
clearvars
%addpath Z:\spin-QED\Theory_matlab\Matlab_functions

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C 

A=117e6;
B0=0.2*5;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+600e6;%-3e8;%6e9; %double-dot tunel rate
Edc=000;%-110;%-330;
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
H_Edc=tensor(identity,e*Edc*L/h*sigma_x,identity);
H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;

[H_eVec,H_eVal]=eig(H_0);

% find correct basis
[Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edc*L/h*sigma_x);
UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);%H_eVec(:,3);%
DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);%H_eVec(:,1);%
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
nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb));

%% find adiabatic pulse shape Eac(t)
clear Eac H_expm
dt=.5e-10;
K=30;
Tcoarse = 200;
PsiDim = 1000;

eo=sqrt(Vt^2+(e*Edc*2*L/h)^2);
DE=(eo-nuE)*pi;

Psit = zeros(8,PsiDim+1);
PsitB = zeros(1,PsiDim+1);
Psit(:,1)=UP;
PsitB(1)=1;

kk=1;
phase=0;
phase1=0;

while phase1<=pi/4
    if kk==1
        Eac(kk)=0;
    else
        Eac(kk)=Eac(kk-1)+dgE/(e*L/2/h*Vt/eo*2*pi);
    end
    gE=e*L/2/h*Vt/eo*Eac(kk)*2*pi;
    dgE=(K^-1)*(DE^2+gE^2)^1.5*dt/DE; %IS THAT TRUE?
    
    % find phase as a function of Eac and time
    H_expm_floquet=eye(8);
    for tt=1:Tcoarse
        H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk)*L/h*sigma_x,identity);
        H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
    end
    H_expm(:,:,kk)=H_expm_floquet^(dt*nuE);
    
    if kk~=1
        tmax=(eo-nuE)/(e*Eac(kk)*2*L*A/h/16*Vt^2/eo^2)*40;
        tdim1=tmax*nuE;
        H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
        for t=1:PsiDim
            Psit(:,t+1)=H_expm_dyn*Psit(:,t);
            PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
        end
        AA=fft(PsitB);
        [pk,loc] = max(AA(2:PsiDim/2));
        phase1=phase+2*pi*dt*loc/tmax
        if phase1<=pi/4
            phase=phase1
        end
    end
    
    kk=kk+1
end

Eac(kk-1)=Eac(kk-2);
H_expm_floquet=eye(8);
for tt=1:Tcoarse
    H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-2)*L/h*sigma_x,identity);
    H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
end
tmax=(eo-nuE)/(e*Eac(kk-2)*2*L*A/h/16*Vt^2/eo^2)*40;
tdim1=tmax*nuE;
H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
for t=1:PsiDim
    Psit(:,t+1)=H_expm_dyn*Psit(:,t);
    PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
end
AA=fft(PsitB);
[pk,loc] = max(AA(2:PsiDim/2));
dtp=(phase1-pi/4)*tmax/2/pi/loc
H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);

time=[0:dt:(kk-2)*dt ((kk-2)*dt+dtp) ((kk-2)*dt+2*dtp) ((kk-2)*dt+2*dtp+dt):dt:((kk-2)*2*dt+2*dtp)];

for k=1:(kk-1)
    Eac(kk-1+k)=Eac(kk-k);
end
Eac(2*kk-1)=0;

figure(3489);
subplot(2,1,1);
plot(time,Eac);

%% time dynamics
charup=tensor(eye(2),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(2));
chardown=tensor(eye(2),Orb_eVec(:,1)*Orb_eVec(:,1)',eye(2));
down=tensor([1 0]'*[1 0],eye(2),[0 1]'*[0 1]);
up=tensor([0 1]'*[0 1],eye(2),[1 0]'*[1 0]);
int=tensor(eye(2),[1 -1]'*[1 -1]./2,eye(2));
don=tensor(eye(2),[1 1]'*[1 1]./2,eye(2));

Psit=zeros(8,1);
Psit(:,1)=(1*UP+0*DOWN)/sqrt(1);

FFqubit=zeros(1,2*(kk-1));
position=zeros(1,2*(kk-1));
ge=zeros(1,2*(kk-1));
ge(1)=Psit(:,1)'*charup*Psit(:,1)-Psit(:,1)'*chardown*Psit(:,1);
FFqubit(1)=Psit(:,1)'*up*Psit(:,1)-Psit(:,1)'*down*Psit(:,1);
position(1)=Psit(:,1)'*don*Psit(:,1)-Psit(:,1)'*int*Psit(:,1);

Operator=eye(8);

for t=1:(kk-1)
    Operator=H_expm(:,:,t)*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*charup*Psit(:,1+t)-Psit(:,1+t)'*chardown*Psit(:,1+t);
    FFqubit(1+t)=Psit(:,1+t)'*up*Psit(:,1+t)-Psit(:,1+t)'*down*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*don*Psit(:,1+t)-Psit(:,1+t)'*int*Psit(:,1+t);
end

for t=kk:(2*kk-2)
    Operator=H_expm(:,:,2*kk-t-1)*Operator;
    Psit(:,1+t)=Operator*Psit(:,1);
    ge(1+t)=Psit(:,1+t)'*charup*Psit(:,1+t)-Psit(:,1+t)'*chardown*Psit(:,1+t);
    FFqubit(1+t)=Psit(:,1+t)'*up*Psit(:,1+t)-Psit(:,1+t)'*down*Psit(:,1+t);
    position(1+t)=Psit(:,1+t)'*don*Psit(:,1+t)-Psit(:,1+t)'*int*Psit(:,1+t);
end

% fidelities
Basis(:,1)=UP;
Basis(:,2)=DOWN;
Basis(:,3)=(UP+DOWN)/sqrt(2);
Basis(:,4)=(1i*UP+DOWN)/sqrt(2);

phase=angle(Basis(:,1)'*Operator*Basis(:,1))-angle(Basis(:,2)'*Operator*Basis(:,2));
OperId=expm(1i*phase/2*sigma_z)*expm(-1i*pi/4*sigma_x);
%OperId=expm(-1i*pi*tpi2*nuE/2*sigma_z)*expm(-1i*pi/4*sigma_x);%[exp(-1i*pi*tpi2*nuE) 1; 1 exp(+1i*pi*tpi2*nuE)]/sqrt(2);%OperId64=zeros(64);
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
Fid4

figure(3489)
subplot(2,1,2);
plot(time*1e9,FFqubit,time*1e9,position,time*1e9,ge);
legend('FFqubit','position','Charge')

%%
FigHandle = figure(3489);
set(FigHandle, 'Position', [500, 100, 900, 900]);

subplot(2,1,1);
plot(time*1e9,Eac,'LineWidth',2);
ylabel('Eac (V/m)')
xlim([-.3 23.8])
ylim([-10 410])
set(gca,'XTick',[0 10 20],'YTick',[0 200 400],'Position', [0.15 0.65 0.8 0.3],'FontSize',50,'XTickLabel','','ticklength',[0.02 0.1])

subplot(2,1,2);
plot(time*1e9,FFqubit,time*1e9,position,time*1e9,ge);
legend('FFqubit','position','Charge')
xlabel('Time (ns)')
xlim([-.3 23.8])
ylim([-1.1 1.1])
set(gca,'XTick',[0 10 20],'YTick',[-1 0 1],'Position', [0.15 0.15 0.8 0.45],'FontSize',50,'ticklength',[0.02 0.1])

%% sweep Edc noise and average        
Edcnoisedim=3;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!
Edcnoise=-1:2/(Edcnoisedim-1):1;
RMSnoise=100;
Edcnoisemax=RMSnoise/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;
Fid4 = zeros(1,Edcnoisedim);

for ll=1:Edcnoisedim
    H_Edc=tensor(identity,e*Edcnoise(ll)*L/h*sigma_x,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
    [H_eVec,H_eVal]=eig(H_0);
    
    % find correct basis
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z+e*Edcnoise(ll)*L/h*sigma_x);
    UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);%H_eVec(:,3);%
    DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);%H_eVec(:,1);%
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

    Operator=eye(8);
    
    for t=1:(kk-2)
        H_expm_floquet=eye(8);
        for tt=1:Tcoarse
            H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(t)*L/h*sigma_x,identity);
            H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
        end
        H_expm(:,:,t)=H_expm_floquet^(dt*nuE);
        Operator=H_expm(:,:,t)*Operator;
    end
    
    H_expm_floquet=eye(8);
    for tt=1:Tcoarse
        H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-1)*L/h*sigma_x,identity);
        H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
    end
    H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
    Operator=H_expm(:,:,kk-1)*Operator;
    
    for t=kk:(2*kk-2)
        Operator=H_expm(:,:,2*kk-t-1)*Operator;
    end
    
    % fidelities
    Basis(:,1)=UP;
    Basis(:,2)=DOWN;
    Basis(:,3)=(UP+DOWN)/sqrt(2);
    Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
    
    %OperId=expm(1i*phase/2*sigma_z)*expm(-1i*pi/4*sigma_x);
    OperId8=zeros(8,8);
    for ii=1:2
        for jj=1:2
            OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
        end
    end
    
    Fid4(ll) = 0;
    for ii=1:4
        Fid4(ll)=Fid4(ll)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
    end
    Fid4(ll)
end
FidNoiseAv=sum(Fid4)/Edcnoisedim

figure(2345)
plot(Edcnoise,Fid4)


%% Gate fidelity with noise  - Vt vs Edc space
profile on
clear vars;

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C 

A=117e6;
B0=0.2*5;
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
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

Vtdim=18;
detmin=50e6;
detmax=900e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

Edcnoisedim=5;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=100/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

EdcDim=21;
EdcMin=-2000;
EdcMax=2000;
Edc=EdcMin:(EdcMax-EdcMin)/(EdcDim-1):EdcMax;

dt=0.5e-10;
K=30;
Tcoarse = 200;
PsiDim = 1000;

Psit = zeros(8,PsiDim+1);
PsitB = zeros(1,PsiDim+1);
PsitB(1)=1;

tpi2 = zeros(Vtdim,EdcDim);
Fid4 = zeros(Vtdim,EdcDim);
Fid4av = zeros(Vtdim,EdcDim);

for ll = 1:Vtdim
    H_tunel=tensor(identity,Vt(ll)/2*sigma_z,identity);
    for mm = 1:EdcDim
        tic;
        H_Edc=tensor(identity,e*Edc(mm)*L/h*sigma_x,identity);
        H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
        
        %find eigenvectors and nuE
        [H_eVec,H_eVal]=eig(H_0);
        [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc(mm)*L/h*sigma_x);
        UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);%H_eVec(:,3);%
        DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);%H_eVec(:,1);%
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
        nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
        Psit(:,1)=UP;
        
        % find adiabatic pulse shape Eac(t)
        clear Eac H_expm
        eo=sqrt(Vt(ll)^2+(e*Edc(mm)*2*L/h)^2);
        DE=(eo-nuE)*pi;
                
        kk=1;
        phase=0;
        phase1=0;
        
        while phase1<=pi/4
            if kk==1
                Eac(kk)=0;
            else
                Eac(kk)=Eac(kk-1)+dgE/(e*L/2/h*Vt(ll)/eo*2*pi);
            end
            gE=e*L/2/h*Vt(ll)/eo*Eac(kk)*2*pi;
            dgE=(K^-1)*(DE^2+gE^2)^1.5*dt/DE; %IS THAT TRUE?
            
            % find phase as a function of Eac and time
            H_expm_floquet=eye(8);
            for tt=1:Tcoarse
                H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk)*L/h*sigma_x,identity);
                H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
            end
            H_expm(:,:,kk)=H_expm_floquet^(dt*nuE);
            
            if kk~=1
                tmax=(eo-nuE)/(e*Eac(kk)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
                tdim1=tmax*nuE;
                H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
                for t=1:PsiDim
                    Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                    PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
                end
                AA=fft(PsitB);
                [pk,loc] = max(AA(2:PsiDim/2));
                phase1=phase+2*pi*dt*loc/tmax;
                if phase1<=pi/4
                    phase=phase1;
                end
            end
            kk=kk+1;
        end
        
        Eac(kk-1)=Eac(kk-2);
        H_expm_floquet=eye(8);
        for tt=1:Tcoarse
            H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-2)*L/h*sigma_x,identity);
            H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
        end
        tmax=(eo-nuE)/(e*Eac(kk-2)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
        tdim1=tmax*nuE;
        H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
        for t=1:PsiDim
            Psit(:,t+1)=H_expm_dyn*Psit(:,t);
            PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
        end
        AA=fft(PsitB);
        [pk,loc] = max(AA(2:PsiDim/2));
        dtp=(phase1-pi/4)*tmax/2/pi/loc;
        H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);

        tpi2(ll,mm)=(kk-2)*2*dt+2*dtp
        
        % time dynamics
        Operator=eye(8);
        for t=1:(kk-1)
            Operator=H_expm(:,:,t)*Operator;
        end
        
        for t=kk:(2*kk-2)
            Operator=H_expm(:,:,2*kk-t-1)*Operator;
        end

        % fidelities
        Basis(:,1)=UP;
        Basis(:,2)=DOWN;
        Basis(:,3)=(UP+DOWN)/sqrt(2);
        Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
        phase=angle(Basis(:,1)'*Operator*Basis(:,1))-angle(Basis(:,2)'*Operator*Basis(:,2));
        OperId=expm(1i*phase/2*sigma_z)*expm(-1i*pi/4*sigma_x);
        OperId8=zeros(8,8);
        for ii=1:2
            for jj=1:2
                OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
            end
        end
        Fid4(ll,mm) = 0;
        for ii=1:4
            Fid4(ll,mm)=Fid4(ll,mm)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
        end
        Fid4(ll,mm)=Fid4(ll,mm)
        
        % sweep Edc noise and average
        Fid4noise = zeros(1,Edcnoisedim);
        Fid4noise(1,(Edcnoisedim-1)/2+1)=Fid4(ll,mm);
        for nnn=[1:(Edcnoisedim-1)/2 ((Edcnoisedim+1)/2+1):Edcnoisedim]
            H_Edc=tensor(identity,e*(Edc(mm)+Edcnoise(nnn))*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
            
            %find correct basis
            [H_eVec,H_eVal]=eig(H_0);
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*(Edc(mm)+Edcnoise(nnn))*L/h*sigma_x);
            UP=tensor([0 1],Orb_eVec(:,1)',[1 0]);%H_eVec(:,3);%
            DOWN=tensor([1 0],Orb_eVec(:,1)',[0 1]);%H_eVec(:,1);%
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
            % time dynamics
            Operator=eye(8);
            for t=1:(kk-2)
                H_expm_floquet=eye(8);
                for tt=1:Tcoarse
                    H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(t)*L/h*sigma_x,identity);
                    H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
                end
                H_expm(:,:,t)=H_expm_floquet^(dt*nuE);
                Operator=H_expm(:,:,t)*Operator;
            end
            
            H_expm_floquet=eye(8);
            for tt=1:Tcoarse
                H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-1)*L/h*sigma_x,identity);
                H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
            end
            H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
            Operator=H_expm(:,:,kk-1)*Operator;
            
            for t=kk:(2*kk-2)
                Operator=H_expm(:,:,2*kk-t-1)*Operator;
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
            Fid4noise(nnn) = 0;
            for ii=1:4
                Fid4noise(nnn)=Fid4noise(nnn)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
            end
        end
        Fid4av(ll,mm)=(sum(Fid4noise))/Edcnoisedim
        toc;
        save 1Q-fid_adiabatic_04.mat
    end
end
profile viewer

%% Origin
VVVt=Vt'/1e9;
EEdc=-Edc;
TotE=1-Fid4av;
min(min(TotE))
max(max(TotE))

%%
figure(672);
subplot(1,3,1);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,tpi2*1e9);axis xy;colorbar;
title('Gate time (ns)')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(1,3,2);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid4))/log(10));axis xy;colorbar;
title('log_{10}(Bare gate error)')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(1,3,3);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-Fid4av)/log(10));axis xy;colorbar;
title('log_{10}(Gate error with noise)')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')


%% Gate fidelity with noise  - Vt space
profile on
clearvars

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C

A=117e6;
B0=0.2*5;
L=7.5e-9; %half interdot distance
dg=-0.002;%relative change in g-factor (delta_g/g)
Edc=-50;

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
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);
H_Edc=tensor(identity,e*Edc*L/h*sigma_x,identity);

Vtdim=18;
detmin=50e6;
detmax=900e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

Edcnoisedim=5;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=100/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

dt=.5e-10;
K=30;
Tcoarse = 200;
PsiDim = 1000;

Psit = zeros(8,PsiDim+1);
PsitB = zeros(1,PsiDim+1);
PsitB(1)=1;

tpi2 = zeros(1,Vtdim);
Fid4 = zeros(1,Vtdim);
Fid4av = zeros(1,Vtdim);
Power = zeros(1,Vtdim);

for ll = 1:Vtdim
    tic;
    H_tunel=tensor(identity,Vt(ll)/2*sigma_z,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
    
    %find eigenvectors and nuE
    [H_eVec,H_eVal]=eig(H_0);
    [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc*L/h*sigma_x);
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
    nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
    Psit(:,1)=UP;
    
    % find adiabatic pulse shape Eac(t)
    clear Eac H_expm
    eo=sqrt(Vt(ll)^2+(e*Edc*2*L/h)^2);
    DE=(eo-nuE)*pi;
    
    kk=1;
    phase=0;
    phase1=0;
    
    while phase1<=pi/4
        if kk==1
            Eac(kk)=0;
        else
            Eac(kk)=Eac(kk-1)+dgE/(e*L/2/h*Vt(ll)/eo*2*pi);
        end
        gE=e*L/2/h*Vt(ll)/eo*Eac(kk)*2*pi;
        dgE=(K^-1)*(DE^2+gE^2)^1.5*dt/DE; %IS THAT TRUE?
        
        % find phase as a function of Eac and time
        H_expm_floquet=eye(8);
        for tt=1:Tcoarse
            H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk)*L/h*sigma_x,identity);
            H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
        end
        H_expm(:,:,kk)=H_expm_floquet^(dt*nuE);
        
        if kk~=1
            tmax=(eo-nuE)/(e*Eac(kk)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
            tdim1=tmax*nuE;
            H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
            for t=1:PsiDim
                Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
            end
            AA=fft(PsitB);
            [pk,loc] = max(AA(2:PsiDim/2));
            phase1=phase+2*pi*dt*loc/tmax;
            if phase1<=pi/4
                phase=phase1;
            end
        end
        kk=kk+1;
    end
    
    Eac(kk-1)=Eac(kk-2);
    H_expm_floquet=eye(8);
    for tt=1:Tcoarse
        H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-2)*L/h*sigma_x,identity);
        H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
    end
    tmax=(eo-nuE)/(e*Eac(kk-2)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
    tdim1=tmax*nuE;
    H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
    for t=1:PsiDim
        Psit(:,t+1)=H_expm_dyn*Psit(:,t);
        PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
    end
    AA=fft(PsitB);
    [pk,loc] = max(AA(2:PsiDim/2));
    dtp=(phase1-pi/4)*tmax/2/pi/loc;
    H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
    
    tpi2(ll)=(kk-2)*2*dt+2*dtp
    Power(ll)=(sum((Eac(1:kk-2)/1e7).^2)/(kk-2+dtp/dt)+(Eac(kk-1)/1e7)^2*dtp/((kk-2)*dt+dtp))/50
    
    % time dynamics
    Operator=eye(8);
    for t=1:(kk-1)
        Operator=H_expm(:,:,t)*Operator;
    end
    
    for t=kk:(2*kk-2)
        Operator=H_expm(:,:,2*kk-t-1)*Operator;
    end
    
    % fidelities
    Basis(:,1)=UP;
    Basis(:,2)=DOWN;
    Basis(:,3)=(UP+DOWN)/sqrt(2);
    Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
    phase=angle(Basis(:,1)'*Operator*Basis(:,1))-angle(Basis(:,2)'*Operator*Basis(:,2));
    OperId=expm(1i*phase/2*sigma_z)*expm(-1i*pi/4*sigma_x);
    OperId8=zeros(8,8);
    for ii=1:2
        for jj=1:2
            OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
        end
    end
    Fid4(ll) = 0;
    for ii=1:4
        Fid4(ll)=Fid4(ll)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
    end
    Fid4(ll)=Fid4(ll)
    
    % sweep Edc noise and average
    Fid4noise = zeros(1,Edcnoisedim);
    Fid4noise(1,(Edcnoisedim-1)/2+1)=Fid4(ll);
    for nnn=[1:(Edcnoisedim-1)/2 ((Edcnoisedim+1)/2+1):Edcnoisedim]
        H_Edc=tensor(identity,e*(Edc+Edcnoise(nnn))*L/h*sigma_x,identity);
        H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
        
        %find correct basis
        [H_eVec,H_eVal]=eig(H_0);
        [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*(Edc+Edcnoise(nnn))*L/h*sigma_x);
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
        % time dynamics
        Operator=eye(8);
        for t=1:(kk-2)
            H_expm_floquet=eye(8);
            for tt=1:Tcoarse
                H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(t)*L/h*sigma_x,identity);
                H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
            end
            H_expm(:,:,t)=H_expm_floquet^(dt*nuE);
            Operator=H_expm(:,:,t)*Operator;
        end
        
        H_expm_floquet=eye(8);
        for tt=1:Tcoarse
            H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-1)*L/h*sigma_x,identity);
            H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
        end
        H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
        Operator=H_expm(:,:,kk-1)*Operator;
        
        for t=kk:(2*kk-2)
            Operator=H_expm(:,:,2*kk-t-1)*Operator;
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
        Fid4noise(nnn) = 0;
        for ii=1:4
            Fid4noise(nnn)=Fid4noise(nnn)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
        end
    end
    Fid4av(ll)=(sum(Fid4noise))/Edcnoisedim
    toc;
    save 1Q-fid_adiabatic_Vt_03.mat
end
profile viewer

%% Origin
VVt=Vt'/1e9;
Gtime=tpi2'*1e9;
LeakN=1-Fid4';
TotN=1-Fid4av';
DPower=Power';

%%
figure(67381);
subplot(2,2,1);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,tpi2*1e9);
ylabel('Gate time (ns)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,2);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid4))/log(10));
ylabel('log_{10}(Bare gate error)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,3);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-Fid4av)/log(10));
ylabel('log_{10}(Gate error with noise)')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,4);
plot((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,Power);
ylabel('Power (W)')
xlabel('V_t - \epsilon_{ff} (MHz)')
 
%%
figure(67371);
subplot(1,2,1);
AX = plotyy((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(Power)/log(10),(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,tpi2*1e9);%log()/log(10);
set(AX(2),'yscale','log')
set(AX(2),'ytickmode','auto')
AX(2).XTickMode = 'auto'
AX(2).XAxisLocation = 'top'
AX(1).Box = 'off'
ylabel('drive power (W), gate time (ns)')
xlabel('\delta_{so} (MHz)')
subplot(1,2,2);
AX = plotyy((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,1-Fid4,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,1-Fid4av);%log()/log(10);
set(AX(2),'yscale','log')
set(AX(2),'ytickmode','auto')
AX(2).XTickMode = 'auto'
AX(2).XAxisLocation = 'top'
AX(1).Box = 'off'
ylabel('gate error')
xlabel('\delta_{so} (MHz)')


%% Gate fidelity with noise  - Vt vs K space
profile on
clearvars

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C

A=117e6;
B0=0.2*5;
L=7.5e-9; %half interdot distance
dg=-0.002;%relative change in g-factor (delta_g/g)
Edc=-50;

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
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

Vtdim=6;
detmin=300e6;
detmax=900e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

RMSnoise=1000;
Edcnoisedim=3;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=RMSnoise/sqrt(sum(Edcnoise.^2)/Edcnoisedim);
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

Kdim=6;
Kmin=10;
Kmax=100;
K=logspace(log(Kmin)/log(10),log(Kmax)/log(10),Kdim);

dt=.5e-10;
Tcoarse = 200;
PsiDim = 1000;

Psit = zeros(8,PsiDim+1);
PsitB = zeros(1,PsiDim+1);
PsitB(1)=1;

tpi2 = zeros(Kdim,Vtdim);
Fid4 = zeros(Kdim,Vtdim);
Fid4av = zeros(Kdim,Vtdim);
Power = zeros(Kdim,Vtdim);

for ll = 1:Vtdim
    tic;
    H_tunel=tensor(identity,Vt(ll)/2*sigma_z,identity);
    
    for kkkk=1:Kdim
        %find eigenvectors and nuE
        H_Edc=tensor(identity,e*Edc*L/h*sigma_x,identity);
        H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
        [H_eVec,H_eVal]=eig(H_0);
        [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc*L/h*sigma_x);
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
        nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
        Psit(:,1)=UP;
        
        % find adiabatic pulse shape Eac(t)
        clear Eac H_expm
        eo=sqrt(Vt(ll)^2+(e*Edc*2*L/h)^2);
        DE=(eo-nuE)*pi;
        
        kk=1;
        phase=0;
        phase1=0;
        
        while phase1<=pi/4
            if kk==1
                Eac(kk)=0;
            else
                Eac(kk)=Eac(kk-1)+dgE/(e*L/2/h*Vt(ll)/eo*2*pi);
            end
            gE=e*L/2/h*Vt(ll)/eo*Eac(kk)*2*pi;
            dgE=(K(kkkk)^-1)*(DE^2+gE^2)^1.5*dt/DE; %IS THAT TRUE?
            
            % find phase as a function of Eac and time
            H_expm_floquet=eye(8);
            for tt=1:Tcoarse
                H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk)*L/h*sigma_x,identity);
                H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
            end
            H_expm(:,:,kk)=H_expm_floquet^(dt*nuE);
            
            if kk~=1
                tmax=(eo-nuE)/(e*Eac(kk)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
                tdim1=tmax*nuE;
                H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
                for t=1:PsiDim
                    Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                    PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
                end
                AA=fft(PsitB);
                [pk,loc] = max(AA(2:PsiDim/2));
                phase1=phase+2*pi*dt*loc/tmax;
                if phase1<=pi/4
                    phase=phase1;
                end
            end
            kk=kk+1;
        end
        
        Eac(kk-1)=Eac(kk-2);
        H_expm_floquet=eye(8);
        for tt=1:Tcoarse
            H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-2)*L/h*sigma_x,identity);
            H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
        end
        tmax=(eo-nuE)/(e*Eac(kk-2)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
        tdim1=tmax*nuE;
        H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
        for t=1:PsiDim
            Psit(:,t+1)=H_expm_dyn*Psit(:,t);
            PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
        end
        AA=fft(PsitB);
        [pk,loc] = max(AA(2:PsiDim/2));
        dtp=(phase1-pi/4)*tmax/2/pi/loc;
        H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
        
        tpi2(kkkk,ll)=(kk-2)*2*dt+2*dtp
        Power(kkkk,ll)=(sum((Eac(1:kk-2)/1e7).^2)/(kk-2+dtp/dt)+(Eac(kk-1)/1e7)^2*dtp/((kk-2)*dt+dtp))/50
        
        % time dynamics
        Operator=eye(8);
        for t=1:(kk-1)
            Operator=H_expm(:,:,t)*Operator;
        end
        
        for t=kk:(2*kk-2)
            Operator=H_expm(:,:,2*kk-t-1)*Operator;
        end
        
        % fidelities
        Basis(:,1)=UP;
        Basis(:,2)=DOWN;
        Basis(:,3)=(UP+DOWN)/sqrt(2);
        Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
        phase=angle(Basis(:,1)'*Operator*Basis(:,1))-angle(Basis(:,2)'*Operator*Basis(:,2));
        OperId=expm(1i*phase/2*sigma_z)*expm(-1i*pi/4*sigma_x);
        OperId8=zeros(8,8);
        for ii=1:2
            for jj=1:2
                OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
            end
        end
        Fid4(kkkk,ll) = 0;
        for ii=1:4
            Fid4(kkkk,ll)=Fid4(kkkk,ll)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
        end
        Fid4(kkkk,ll)=Fid4(kkkk,ll)
        
        % sweep Edc noise and average
        Fid4noise = zeros(1,Edcnoisedim);
        Fid4noise((Edcnoisedim-1)/2+1)=Fid4(kkkk,ll);
        for nnn=[1:(Edcnoisedim-1)/2 ((Edcnoisedim+1)/2+1):Edcnoisedim]
            H_Edc=tensor(identity,e*(Edc+Edcnoise(nnn))*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
            
            %find correct basis
            [H_eVec,H_eVal]=eig(H_0);
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*(Edc+Edcnoise(nnn))*L/h*sigma_x);
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
            % time dynamics
            Operator=eye(8);
            for t=1:(kk-2)
                H_expm_floquet=eye(8);
                for tt=1:Tcoarse
                    H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(t)*L/h*sigma_x,identity);
                    H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
                end
                H_expm(:,:,t)=H_expm_floquet^(dt*nuE);
                Operator=H_expm(:,:,t)*Operator;
            end
            
            H_expm_floquet=eye(8);
            for tt=1:Tcoarse
                H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-1)*L/h*sigma_x,identity);
                H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
            end
            H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
            Operator=H_expm(:,:,kk-1)*Operator;
            
            for t=kk:(2*kk-2)
                Operator=H_expm(:,:,2*kk-t-1)*Operator;
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
            Fid4noise(nnn) = 0;
            for ii=1:4
                Fid4noise(nnn)=Fid4noise(nnn)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
            end
        end
        Fid4av(kkkk,ll)=(sum(Fid4noise))/Edcnoisedim
        toc;
        save 1Q-fid_adiabatic_Vt_K_03_1000VmRMS.mat
    end
end
profile viewer

%%
figure(67253);
subplot(2,2,1);
imagesc((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(K)/log(10),tpi2*1e9);axis xy;colorbar;
title('Gate time (ns)')
ylabel('log_{10}K')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,2);
imagesc((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(K)/log(10),log(abs(1-Fid4))/log(10));axis xy;colorbar;
title('log_{10}(Bare gate error)')
ylabel('log_{10}K')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,3);
imagesc((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(K)/log(10),log(1-Fid4av)/log(10));axis xy;colorbar;
title('log_{10}(Gate error with noise)')
ylabel('log_{10}K')
xlabel('V_t - \epsilon_{ff} (MHz)')
subplot(2,2,4);
imagesc((Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(K)/log(10),log(Power)/log(10));axis xy;colorbar;
title('log_{10}(Power) (W)')
ylabel('log_{10}K')
xlabel('V_t - \epsilon_{ff} (MHz)')


%% Gate fidelity with noise  - sweep Vt, K, Edc and RMSnoise - K vs RMSnoise space
profile on
clearvars

h=6.62606957e-34; %J*s
e=1.602176565e-19; %C

A=117e6;
B0=0.2*5;
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
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

Edcdim=5;
Edcmin=-1000;
Edcmax=300;
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;

Vtdim=20;
detmin=50e6;
detmax=2500e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(Vtdim-1):detmax); %double-dot tunel rate

RMSdim=9;
RMSmin=10;
RMSmax=1000;
RMSnoise=logspace(log(RMSmin)/log(10),log(RMSmax)/log(10),RMSdim);

Edcnoisedim=5;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!

Kdim=9;
Kmin=3;
Kmax=300;
K=logspace(log(Kmin)/log(10),log(Kmax)/log(10),Kdim);

dtrough=.5e-9;
Tcoarse = 200;
PsiDim = 1000;

Psit = zeros(8,PsiDim+1);
PsitB = zeros(1,PsiDim+1);
PsitB(1)=1;

tpi2 = zeros(Edcdim,Vtdim,Kdim);
Fid4 = zeros(Edcdim,Vtdim,Kdim,RMSdim,Edcnoisedim);

for ll = 1:Vtdim
    H_tunel=tensor(identity,Vt(ll)/2*sigma_z,identity);
    for kkkk=1:Kdim
        for eed=1:Edcdim
            tic;
            %find eigenvectors and nuE
            H_Edc=tensor(identity,e*Edc(eed)*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
            [H_eVec,H_eVal]=eig(H_0);
            [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*Edc(eed)*L/h*sigma_x);
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
            nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
            Psit(:,1)=UP;
            
            % find approximate tpi2
            clear Eac H_expm
            eo=sqrt(Vt(ll)^2+(e*Edc(eed)*2*L/h)^2);
            DE=(eo-nuE)*pi;
            
            kk=1;
            phase=0;
            
            while phase<=pi/4
                if kk==1
                    Eac(kk)=0;
                else
                    Eac(kk)=Eac(kk-1)+dgE/(e*L/2/h*Vt(ll)/eo*2*pi);
                end
                gE=e*L/2/h*Vt(ll)/eo*Eac(kk)*2*pi;
                dgE=(K(kkkk)^-1)*(DE^2+gE^2)^1.5*dtrough/DE; %IS THAT TRUE?
                
                % find phase as a function of Eac and time
                H_expm_floquet=eye(8);
                for tt=1:Tcoarse
                    H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk)*L/h*sigma_x,identity);
                    H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
                end
                H_expm(:,:,kk)=H_expm_floquet^(dtrough*nuE);
                
                if kk~=1
                    tmax=(eo-nuE)/(e*Eac(kk)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
                    tdim1=tmax*nuE;
                    H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
                    for t=1:PsiDim
                        Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                        PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
                    end
                    AA=fft(PsitB);
                    [pk,loc] = max(AA(2:PsiDim/2));
                    phase=phase+2*pi*dtrough*loc/tmax;
                end
                kk=kk+1;
            end
            
            (kk-1)*2*dtrough
            dt=(kk-1)*dtrough/200;         
            
            % find adiabatic pulse shape Eac(t)
            clear Eac H_expm
            
            kk=1;
            phase=0;
            phase1=0;
            
            while phase1<=pi/4
                if kk==1
                    Eac(kk)=0;
                else
                    Eac(kk)=Eac(kk-1)+dgE/(e*L/2/h*Vt(ll)/eo*2*pi);
                end
                gE=e*L/2/h*Vt(ll)/eo*Eac(kk)*2*pi;
                dgE=(K(kkkk)^-1)*(DE^2+gE^2)^1.5*dt/DE; %IS THAT TRUE?
                
                % find phase as a function of Eac and time
                H_expm_floquet=eye(8);
                for tt=1:Tcoarse
                    H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk)*L/h*sigma_x,identity);
                    H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
                end
                H_expm(:,:,kk)=H_expm_floquet^(dt*nuE);
                
                if kk~=1
                    tmax=(eo-nuE)/(e*Eac(kk)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
                    tdim1=tmax*nuE;
                    H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
                    for t=1:PsiDim
                        Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                        PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
                    end
                    AA=fft(PsitB);
                    [pk,loc] = max(AA(2:PsiDim/2));
                    phase1=phase+2*pi*dt*loc/tmax;
                    if phase1<=pi/4
                        phase=phase1;
                    end
                end
                kk=kk+1;
            end
            
            Eac(kk-1)=Eac(kk-2);
            H_expm_floquet=eye(8);
            for tt=1:Tcoarse
                H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-2)*L/h*sigma_x,identity);
                H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
            end
            tmax=(eo-nuE)/(e*Eac(kk-2)*2*L*A/h/16*Vt(ll)^2/eo^2)*40;
            tdim1=tmax*nuE;
            H_expm_dyn=H_expm_floquet^(tdim1/PsiDim);
            for t=1:PsiDim
                Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                PsitB(t+1)=abs(Psit(:,t+1)'*UP)^2;
            end
            AA=fft(PsitB);
            [pk,loc] = max(AA(2:PsiDim/2));
            dtp=(phase1-pi/4)*tmax/2/pi/loc;
            H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
            
            tpi2(eed,ll,kkkk)=(kk-2)*2*dt+2*dtp;
            tpi2(eed,ll,kkkk)
            
            % time dynamics
            Operator=eye(8);
            for t=1:(kk-1)
                Operator=H_expm(:,:,t)*Operator;
            end
            
            for t=kk:(2*kk-2)
                Operator=H_expm(:,:,2*kk-t-1)*Operator;
            end
            
            % fidelities
            Basis(:,1)=UP;
            Basis(:,2)=DOWN;
            Basis(:,3)=(UP+DOWN)/sqrt(2);
            Basis(:,4)=(1i*UP+DOWN)/sqrt(2);
            phase=angle(Basis(:,1)'*Operator*Basis(:,1))-angle(Basis(:,2)'*Operator*Basis(:,2));
            OperId=expm(1i*phase/2*sigma_z)*expm(-1i*pi/4*sigma_x);
            OperId8=zeros(8,8);
            for ii=1:2
                for jj=1:2
                    OperId8=OperId8+OperId(ii,jj)*Basis(:,ii)*Basis(:,jj)';
                end
            end
            %Fid4(eed,:,(Edcnoisedim-1)/2+1,kkkk,ll) = 0;
            for ii=1:4
                Fid4(eed,ll,kkkk,:,(Edcnoisedim+1)/2)=Fid4(eed,ll,kkkk,:,(Edcnoisedim+1)/2)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
            end
            
            % sweep Edc noise and average
            for rms=1:RMSdim
                Edcnoise=-1:2/(Edcnoisedim-1):1;
                Edcnoisemax=RMSnoise(rms)/sqrt(sum(Edcnoise.^2)/Edcnoisedim);
                Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;
                for nnn=[1:(Edcnoisedim-1)/2 ((Edcnoisedim+1)/2+1):Edcnoisedim]
                    H_Edc=tensor(identity,e*(Edc(eed)+Edcnoise(nnn))*L/h*sigma_x,identity);
                    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
                    
                    %find correct basis
                    [H_eVec,H_eVal]=eig(H_0);
                    [Orb_eVec,Orb_eVal]=eig(Vt(ll)/2*sigma_z+e*(Edc(eed)+Edcnoise(nnn))*L/h*sigma_x);
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
                    % time dynamics
                    Operator=eye(8);
                    for t=1:(kk-2)
                        H_expm_floquet=eye(8);
                        for tt=1:Tcoarse
                            H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(t)*L/h*sigma_x,identity);
                            H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
                        end
                        H_expm(:,:,t)=H_expm_floquet^(dt*nuE);
                        Operator=H_expm(:,:,t)*Operator;
                    end
                    
                    H_expm_floquet=eye(8);
                    for tt=1:Tcoarse
                        H=H_0+cos(2*pi*tt/Tcoarse)*tensor(identity,e*Eac(kk-1)*L/h*sigma_x,identity);
                        H_expm_floquet=expm(-1i*2*pi*H/(nuE*Tcoarse))*H_expm_floquet;
                    end
                    H_expm(:,:,kk-1)=H_expm_floquet^(dtp*nuE);
                    Operator=H_expm(:,:,kk-1)*Operator;
                    
                    for t=kk:(2*kk-2)
                        Operator=H_expm(:,:,2*kk-t-1)*Operator;
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
                    %Fid4noise(nnn) = 0;
                    for ii=1:4
                        Fid4(eed,ll,kkkk,rms,nnn)=Fid4(eed,ll,kkkk,rms,nnn)+abs((Operator*Basis(:,ii))'*(OperId8*Basis(:,ii)))^2/4;
                    end
                end
            end
            toc;
            ll,kkkk,eed
            save 1Q-fid_adiabatic_5D_06.mat
        end
    end
end
profile viewer

%% plot gate time and error for an all K and RMSnoise values (only works if RMSdim=Kdim=3)
figure(672531);
subplot(3,5,1);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,tpi2(:,:,1)'*1e9);axis xy;colorbar;
title('time (ns), K=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,2);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid4(:,:,1,1,(Edcnoisedim+1)/2)'))/log(10));axis xy;colorbar;
title('log_{10}(leakage error), K=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,3);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,1,1,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=10, RMS=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,4);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,1,2,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=10, RMS=100')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,5);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,1,3,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=10, RMS=1000')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,1+5);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,tpi2(:,:,2)'*1e9);axis xy;colorbar;
title('time (ns), K=45')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,2+5);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid4(:,:,2,1,(Edcnoisedim+1)/2)'))/log(10));axis xy;colorbar;
title('log_{10}(leakage error), K=45')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,3+5);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,2,1,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=45, RMS=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,4+5);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,2,2,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=45, RMS=100')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,5+5);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,2,3,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=45, RMS=1000')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,1+10);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,tpi2(:,:,3)'*1e9);axis xy;colorbar;
title('time (ns), K=200')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,2+10);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(abs(1-Fid4(:,:,3,1,(Edcnoisedim+1)/2)'))/log(10));axis xy;colorbar;
title('log_{10}(leakage error), K=200')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,3+10);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,3,1,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=200, RMS=10')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,4+10);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,3,2,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=200, RMS=100')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')
subplot(3,5,5+10);
imagesc(-Edc,(Vt-sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4))/1e6,log(1-sum(Fid4(:,:,3,3,:),5)'/Edcnoisedim)/log(10));axis xy;colorbar;
title('log_{10}(total error), K=200, RMS=1000')
xlabel('Edc final (V/m)')
ylabel('V_t - \epsilon_{ff} (MHz)')

%% plot optimized total error vs RMSnoise and K - optimze Edc as well
%tpi2 = zeros(Edcdim,Vtdim,Kdim);
%Fid4 = zeros(Edcdim,Vtdim,Kdim,RMSdim,Edcnoisedim);
TotErr=zeros(Kdim,RMSdim);
VVt=zeros(Kdim,RMSdim);
EEdc=zeros(Kdim,RMSdim);
Gtime=zeros(Kdim,RMSdim);

for kkk=1:Kdim
    for rms=1:RMSdim
        [mm,eedc]=min(1-sum(Fid4(:,:,kkk,rms,:),5)/Edcnoisedim);
        [mm,vvt]=min(mm);
        eedc=eedc(vvt);
        EEdc(kkk,rms)=-Edc(eedc);
        VVt(kkk,rms)=Vt(vvt);
        Gtime(kkk,rms)=tpi2(eedc,vvt,kkk);
        TotErr(kkk,rms)=mm;
    end
end
    
figure(672532);
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
title('optimal Edc (V/m)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')

%% Origin
Rnoise=RMSnoise;
Kk=K';
TTotErr=TotErr;
GGtime=sum(Gtime,2)/RMSdim;

%% plot optimized total error vs RMSnoise and K - assume Edc=0 - MANUALLY SELECT EDC=0 INDEX!!!!!!!!!!!!!
%tpi2 = zeros(Edcdim,Vtdim,Kdim);
%Fid4 = zeros(Edcdim,Vtdim,Kdim,RMSdim,Edcnoisedim);
TotErr=zeros(Kdim,RMSdim);
VVt=zeros(Kdim,RMSdim);
EEdc=zeros(Kdim,RMSdim);
Gtime=zeros(Kdim,RMSdim);

eedc=4; % MANUALLY SELECT EDC=0 INDEX!!!!!!!!!!!!!

for kkk=1:Kdim
    for rms=1:RMSdim
        [mm,vvt]=min(1-sum(Fid4(eedc,:,kkk,rms,:),5)/Edcnoisedim);
        EEdc(kkk,rms)=-Edc(eedc);
        VVt(kkk,rms)=Vt(vvt);
        Gtime(kkk,rms)=tpi2(eedc,vvt,kkk);
        TotErr(kkk,rms)=mm;
    end
end
    
figure(6725321);
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
title('optimal Edc (V/m)')
ylabel('log_{10}K')
xlabel('log_{10}(E_{rms}) (V/m)')