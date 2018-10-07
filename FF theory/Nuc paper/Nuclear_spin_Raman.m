%% Noise nucleus rotating frame, resonant charge transition, sweep Edc - LEVEL DIAGRAM IN ROTATING FRAME - drive perturbs charge
clearvars
addpath Z:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s
k_B= 1.38064852e-23; %m2kgs-2K-1

L=7.5e-9; %half interdot distance
A=117e6;
B0=0.2;%
Bac=0.6e-3;%0.477e-3;%
dg=-0.002; %relative difference in g-factor at interface

Eac = 32*0.0001; % Electric drive
tdim=2000;
T_coarse=200;

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% Basis: nucleus x charge x electron
% 01) u,e,u
% 02) u,e,d
% 03) u,g,u
% 04) u,g,d
% 05) d,e,u
% 06) d,e,d
% 07) d,g,u
% 08) d,g,d

Vt=sqrt(((28e9*0.999+17.2e6)*B0)^2+A^2/4)+6e6+0*300e6;%
nuB=28e9*0.999*B0-A/4-0e6-0*42e6;
%rotating frame! nuB
H_Znuc=tensor((-17.2e6*B0-nuB)/2*sigma_z,identity,identity); 
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z)+tensor(identity,identity,(-nuB)/2*sigma_z);
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_ESRnuc=tensor(-17.2e6*Bac/2/2*sigma_x,identity,identity);
H_ESRel=tensor(identity,identity-(identity/2+sigma_x/2)*0.002*0,28e9*Bac/2/2*sigma_x);

nuc_e=tensor([0 1],[0 1],[0 1]); %excited nuclear state
nuc_g=tensor([1 0],[0 1],[0 1]); %ground nuclear state
ele_e=tensor([0 1],[0 1],[1 0]); %excited electron and nuclear state
char_e=tensor([1 0],[1 0],[0 1]); %excited charge state

Edcdim=300;
Edcmin=-0.6e4*B0/2;
Edcmax=0.6e4*B0/2;
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;

ggB0=zeros(Edcdim,1);
gnB0=zeros(Edcdim,1);
gcB0=zeros(Edcdim,1);
RelTheor=zeros(Edcdim,1);
ensrot=zeros(Edcdim,1);
hyper=zeros(Edcdim,1);
gEns=zeros(Edcdim,1);
Vtp=zeros(Edcdim,1);
alpha=zeros(Edcdim,1);
beta=zeros(Edcdim,1);
alphaup=zeros(Edcdim,1);
betaup=zeros(Edcdim,1);
alphadown1=zeros(Edcdim,1);
betadown1=zeros(Edcdim,1);
alphadown2=zeros(Edcdim,1);
betadown2=zeros(Edcdim,1);
Rel=zeros(Edcdim,1);
nucfrac=zeros(Edcdim,1);
flipfrac=zeros(Edcdim,1);
charfrac=zeros(Edcdim,1);
nucfractheory=zeros(Edcdim,1);

edressed=zeros(Edcdim,1);
eo=zeros(Edcdim,1);
deltaB=zeros(Edcdim,1);
nuRabi=zeros(Edcdim,1);
Psit=zeros(8,tdim+1);
PsitB=zeros(1,tdim+1);

for ii=1:Edcdim
    tic;
    eo(ii,1)=sqrt(Vt^2+(e*Edc(ii)*L*2/h)^2);
    H_Edc=-tensor(identity,e*Edc(ii)*L/h*sigma_x,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
    [H_eVec,H_eVal]=eig(H_0);
    
    fid1=0;
    fid2=0;
    fid3=0;
    fid4=0;
    aa=0;
    bb=0;
    cc=0;
    dd=0;
    for i=1:8
        AA=abs(nuc_e*H_eVec(:,i))^2;
        BB=abs(nuc_g*H_eVec(:,i))^2;
        CC=abs(ele_e*H_eVec(:,i))^2;
        DD=abs(char_e*H_eVec(:,i))^2;
        if AA>fid1
            fid1=AA;
            aa=i;
        end
        if BB>fid2
            fid2=BB;
            bb=i;
        end
        if CC>fid3
            fid3=CC;
            cc=i;
        end
        if DD>fid4
            fid4=DD;
            dd=i;
        end
    end
    
    ggB0(ii,1)=abs(H_eVal(aa,aa)-H_eVal(bb,bb)); %nuclear flip
    gnB0(ii,1)=abs(H_eVal(cc,cc)-H_eVal(bb,bb)); %flip flop
    gcB0(ii,1)=abs(H_eVal(dd,dd)-H_eVal(bb,bb)); %charge change
    
    Vtp(ii,1) = sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2);
    Ap = A/2*(1-((e*2*L*Edc(ii)/h)/Vtp(ii,1)));
    ge = 28e9*(1+(1+((e*2*L*Edc(ii)/h)/Vtp(ii,1)))*dg/2);
    g_ffo=A/4*Vt/Vtp(ii,1);
    gn = 17.2e6;
    %electron,nucleus energies
    Euu =  ((ge - gn)*B0 + Ap/2)/2;
    Edd =  (-1*(ge - gn)*B0 + Ap/2)/2;
    Edu =  (-1*sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    Eud =  (sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    
    %gnB0(ii,1)=Eud-Edu;
    
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z-e*Edc(ii)*L/h*sigma_x);
    charup=tensor(eye(2),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(2));
    %relaxation rate
    Rel(ii,1)=(H_eVec(:,aa)'*charup*H_eVec(:,aa))*2.37e-24*Vtp(ii,1)*Vt^2;
    %(H_eVec(:,aa)'*charup*H_eVec(:,aa))
    
    nuc_e=tensor([0 1],Orb_eVec(:,1)',[0 1]);
    nuc_g=tensor([1 0],Orb_eVec(:,1)',[0 1]);
    ele_e=tensor([0 1],Orb_eVec(:,1)',[1 0]);
    char_e=tensor([1 0],Orb_eVec(:,2)',[0 1]);
    nucfrac(ii,1)=abs(nuc_e*H_eVec(:,aa))^2;
    flipfrac(ii,1)=abs(nuc_e*H_eVec(:,cc))^2;
    %charfrac(ii,1)=abs(char_e*H_eVec(:,aa))^2;
    
    gB=ge*Bac/4;    
    deltaB(ii,1)=(Eud-Edd)-nuB;
    deltaSO=Vtp(ii,1)-(gn*B0+nuB+Ap/2+deltaB(ii,1));
    gSO=A/4*Vt/Vtp(ii,1);
    %excited flip flop hybrid states (19)
    plus=(deltaSO+sqrt(deltaSO^2+(2*gSO)^2))/2/gSO; %phi
    minus=(deltaSO-sqrt(deltaSO^2+(2*gSO)^2))/2/gSO; %theta
    alpha(ii,1)=minus/sqrt(minus^2+1); %alpha
    beta(ii,1)=plus/sqrt(plus^2+1); %beta
    DSO=deltaSO/2*(sqrt(1+(2*gSO/deltaSO)^2)-1); %Ddrive
    %(deltaB-DSO)/gB
    
    %electron flip with nuc up
    plusup=((deltaB(ii,1)+Ap)+sqrt((deltaB(ii,1)+Ap)^2+(2*gB)^2))/2/gB;
    minusup=((deltaB(ii,1)+Ap)-sqrt((deltaB(ii,1)+Ap)^2+(2*gB)^2))/2/gB;
    alphaup(ii,1)=minusup/sqrt(minusup^2+1); %alpha3
    betaup(ii,1)=plusup/sqrt(plusup^2+1); %beta3
    
    %electron flip with nuc down, lower hybrid (Ap missing??)
    plusdown1=((deltaB(ii,1)-DSO)+sqrt((deltaB(ii,1)-DSO)^2+(2*gB*beta(ii,1))^2))/2/gB/beta(ii,1); 
    minusdown1=((deltaB(ii,1)-DSO)-sqrt((deltaB(ii,1)-DSO)^2+(2*gB*beta(ii,1))^2))/2/gB/beta(ii,1);
    alphadown1(ii,1)=minusdown1/sqrt(minusdown1^2+1); %alpha1
    betadown1(ii,1)=plusdown1/sqrt(plusdown1^2+1); %beta1
   
    %electron flip with nuc down, upper hybrid (Ap missing??)
    plusdown2=((deltaB(ii,1)+deltaSO+DSO)+sqrt((deltaB(ii,1)+deltaSO+DSO)^2+(-2*gB*alpha(ii,1))^2))/2/gB/(-alpha(ii,1));
    minusdown2=((deltaB(ii,1)+deltaSO+DSO)-sqrt((deltaB(ii,1)+deltaSO+DSO)^2+(-2*gB*alpha(ii,1))^2))/2/gB/(-alpha(ii,1));
    alphadown2(ii,1)=minusdown2/sqrt(minusdown2^2+1); %alpha2
    betadown2(ii,1)=plusdown2/sqrt(plusdown2^2+1); %beta2
    %alphadown2.^2+betadown2.^2
    
    %relaxation rate
    RelTheor(ii,1)=(min(abs(betadown1(ii,1)),abs(alphadown1(ii,1)))*alpha(ii,1)+alphadown2(ii,1)*beta(ii,1))^2*2.37e-24*Vtp(ii,1)*Vt^2;
    
    hyper(ii,1) = Edd-Edu + nuB;
    %e_ns (22)
    ensrot(ii,1) = gn*B0 + Ap/2 + (deltaB(ii,1)+Ap)/2*(sqrt(1+(2*gB/(deltaB(ii,1)+Ap))^2)-1) - (deltaB(ii,1)-DSO)/2*(sqrt(1+(2*gB*beta(ii,1)/(deltaB(ii,1)-DSO))^2)-1) - (deltaB(ii,1)+deltaSO+DSO)/2*(sqrt(1+(2*gB*(-alpha(ii,1))/(deltaB(ii,1)+deltaSO+DSO))^2)-1);
    %effrot(ii,1) = gn*B0 + sign(abs(betaup)-abs(alphaup))*sqrt((Ap+deltaB)^2+(2*gB)^2)/2 + sign(abs(betadown)-abs(alphadown))*sqrt(deltaB^2+(2*gB)^2)/2;
    
    gE=e*Eac*2*L/4/h*Vt/Vtp(ii,1);
    
    %nuclear coupling strength
    gEns(ii,1)=-gE*betaup(ii,1)*(min(abs(betadown1(ii,1)),abs(alphadown1(ii,1)))*alpha(ii,1)+alphadown2(ii,1)*beta(ii,1));%
    
    nucfractheory(ii,1)=(max(abs(betadown1(ii,1)),abs(alphadown1(ii,1)))*betadown2(ii,1))^2;
    
    % find gate time
    nuE=ggB0(ii,1);
    
    H_expm=eye(8);
    for tt=1:T_coarse
        H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
        H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
    end
    
    tmax=(A/4*Vt/Vtp(ii,1))/(28e9*Bac/4)/(e*Eac*2*L/h/4*Vt/Vtp(ii,1))*100;
    tdim1=tmax*nuE;
    H_expm1=H_expm^(tdim1/tdim);
    
    Psit(:,1)=H_eVec(:,aa);
    for t=1:tdim
        Psit(:,t+1)=H_expm1*Psit(:,t);
        PsitB(t+1)=abs(Psit(:,t+1)'*nuc_g')^2;
    end
    
    AA=fft(PsitB);
    [pk,loc] = max(abs(AA(2:tdim/2)));
    nuRabi(ii,1)=loc/tmax;
    
    ii
    toc;
    edressed(ii,1)=nuB+deltaB(ii,1)+Ap;
end

%%
figure(213)
plot(Edc/1e3,nuB+deltaB,Edc/1e3,eo);%plot(Edc/1e3,nucfrac,Edc/1e3,nucfractheory)
%%
figure('position',[680 364 356 614])
subplot(2,1,1)
plot(Edc/1e3,hyper,Edc/1e3,gcB0,Edc/1e3,gnB0,Edc/1e3,ggB0,Edc/1e3,ensrot+nuB,Edc/1e3,Vtp);%Edc,effrot,%-nuB
%plot(Edc/1e3,hyper,Edc(1:228)/1e3,ggB0(1:228),Edc(229:400)/1e3,ggB0(229:400),Edc(1:228)/1e3,gnB0(1:228),Edc(229:400)/1e3,gnB0(229:400),Edc(1:228)/1e3,ensrot(1:228)+nuB,Edc(229:400)/1e3,ensrot(229:400)+nuB,Edc/1e3,Vtp);%Edc,effrot,%-nuB
legend('eps_ns+vB','charge change','ff','nuc flip','theory','charge')%'theory',
ylim([5.5701e9, 5.627e9]-0*42e6)%-nuB %ylim([min(min(ggB0),min(gnB0)), max(max(ggB0),max(gnB0))])
xlim([-2.5, 2.5])
subplot(4,1,3)
plot(Edc/1e3,nuRabi/2/Eac,Edc/1e3,gEns/Eac)%
legend('numerics','theory')
ylim([-2e3,1.55e5]);
xlim([-2.5, 2.5])
subplot(4,1,4)
plot(Edc/1e3,Rel,Edc/1e3,RelTheor)%
legend('numerics','theory')
ylim([-2e2, 11e3])
xlim([-2.5, 2.5])

%% figure 3a
figure('position',[680 364 356 614])
subplot(2,1,1)
i=1
%weird numbers are for the color only 
plot(Edc(i)/1e3,ggB0(i),'o','markersize',5,'MarkerEdgeColor','none','MarkerFaceColor',[242-161*nucfrac(i) 199-82*nucfrac(i) 92-59*nucfrac(i)]/255)
hold on
plot(Edc(i)/1e3,gnB0(i),'o','markersize',5,'MarkerEdgeColor','none','MarkerFaceColor',[242-161*flipfrac(i) 199-82*flipfrac(i) 92-59*flipfrac(i)]/255)
% for i=2:500
% plot(Edc(i)/1e3,ggB0(i),'o','markersize',5,'MarkerEdgeColor','none','MarkerFaceColor',[242-161*nucfrac(i) 199-82*nucfrac(i) 92-59*nucfrac(i)]/255) %nuclear flip green
% plot(Edc(i)/1e3,gnB0(i),'o','markersize',5,'MarkerEdgeColor','none','MarkerFaceColor',[242-161*flipfrac(i) 199-82*flipfrac(i) 92-59*flipfrac(i)]/255) %ff yellow
% end
hold off
ylim([5.5701e9, 5.627e9])
xlim([-2.5, 2.5])


%%
figure(7352)
plot(Edc,alphaup,Edc,betaup)%
legend('beta1','alpha1')%

%%
figure(7353)
plot(Edc,alphaup,Edc,betadown1,Edc,betadown2,Edc,alpha,Edc,beta,Edc,betaup,Edc,alphadown1,Edc,alphadown2)%
legend('alpha3','beta1','beta2','alpha','beta','beta3','alpha1','alpha2')%
%%
figure(23434)
plot(Edc,nucfrac,Edc,flipfrac,Edc,charfrac,Edc,nucfrac+flipfrac+charfrac)
legend('nuc', 'ff', 'c', 'tot')

%% Noise nucleus rotating frame, resonant charge transition, sweep Edc and nuB
%fig 4
clearvars
addpath Z:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

L=7.5e-9; %half interdot distance
A=117e6;
B0=0.2;%
Bac=0.6e-3;%0.477e-3;%
dg=-0.002; %relative difference in g-factor at interface

Eac = 32*0.0001; % Electric drive
tdim=2000;
T_coarse=200;

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% Basis: nucleus x charge x electron
% 01) u,e,u
% 02) u,e,d
% 03) u,g,u
% 04) u,g,d
% 05) d,e,u
% 06) d,e,d
% 07) d,g,u
% 08) d,g,d

Vt=sqrt(((28e9*0.999+17.2e6)*B0)^2+A^2/4)+6e6+0*300e6;%
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_ESRnuc=tensor(-17.2e6*Bac/2/2*sigma_x,identity,identity);
H_ESRel=tensor(identity,identity-(identity/2+sigma_x/2)*0.002*0,28e9*Bac/2/2*sigma_x);

nuc_e=tensor([0 1],[0 1],[0 1]);
nuc_g=tensor([1 0],[0 1],[0 1]);
ele_e=tensor([0 1],[0 1],[1 0]);
char_e=tensor([1 0],[1 0],[0 1]);

Edcdim=200;
Edcmin=-0.6e4*B0/2;
Edcmax=0.6e4*B0/2;
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;

ggB0=zeros(Edcdim,1);
gnB0=zeros(Edcdim,1);
gcB0=zeros(Edcdim,1);
ensrot=zeros(Edcdim,1);
hyper=zeros(Edcdim,1);
Vtp=zeros(Edcdim,1);
alpha=zeros(Edcdim,1);
beta=zeros(Edcdim,1);
alphaup=zeros(Edcdim,1);
betaup=zeros(Edcdim,1);
alphadown1=zeros(Edcdim,1);
betadown1=zeros(Edcdim,1);
alphadown2=zeros(Edcdim,1);
betadown2=zeros(Edcdim,1);
nucfrac=zeros(Edcdim,1);
flipfrac=zeros(Edcdim,1);
charfrac=zeros(Edcdim,1);
nucfractheory=zeros(Edcdim,1);

Psit=zeros(8,tdim+1);
PsitB=zeros(1,tdim+1);

nuB=(28e9*0.999*B0-A/4-0e6-0*42e6)+(-60e6:2e6:100e6);

gEns=zeros(Edcdim,size(nuB,2));
Rel=zeros(Edcdim,size(nuB,2));
RelTheor=zeros(Edcdim,size(nuB,2));
nuRabi=zeros(Edcdim,size(nuB,2));

for jj=1:size(nuB,2)
    tic;
    H_Znuc=tensor((-17.2e6*B0-nuB(jj))/2*sigma_z,identity,identity);
    H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z)+tensor(identity,identity,(-nuB(jj))/2*sigma_z);

for ii=1:Edcdim
    H_Edc=-tensor(identity,e*Edc(ii)*L/h*sigma_x,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
    [H_eVec,H_eVal]=eig(H_0);
    
    fid1=0;
    fid2=0;
    fid3=0;
    fid4=0;
    aa=0;
    bb=0;
    cc=0;
    dd=0;
    for i=1:8
        AA=abs(nuc_e*H_eVec(:,i))^2;
        BB=abs(nuc_g*H_eVec(:,i))^2;
        CC=abs(ele_e*H_eVec(:,i))^2;
        DD=abs(char_e*H_eVec(:,i))^2;
        if AA>fid1
            fid1=AA;
            aa=i;
        end
        if BB>fid2
            fid2=BB;
            bb=i;
        end
        if CC>fid3
            fid3=CC;
            cc=i;
        end
        if DD>fid4
            fid4=DD;
            dd=i;
        end
    end
    
    ggB0(ii,1)=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
    gnB0(ii,1)=abs(H_eVal(cc,cc)-H_eVal(bb,bb));
    gcB0(ii,1)=abs(H_eVal(dd,dd)-H_eVal(bb,bb));
    
    Vtp(ii,1) = sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2);
    Ap = A/2*(1-((e*2*L*Edc(ii)/h)/Vtp(ii,1)));
    ge = 28e9*(1+(1+((e*2*L*Edc(ii)/h)/Vtp(ii,1)))*dg/2);
    g_ffo=A/4*Vt/Vtp(ii,1);
    gn = 17.2e6;
    %electron,nucleus
    Euu =  ((ge - gn)*B0 + Ap/2)/2;
    Edd =  (-1*(ge - gn)*B0 + Ap/2)/2;
    Edu =  (-1*sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    Eud =  (sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    
    %gnB0(ii,1)=Eud-Edu;
    
    [Orb_eVec,Orb_eVal]=eig(Vt/2*sigma_z-e*Edc(ii)*L/h*sigma_x);
    charup=tensor(eye(2),Orb_eVec(:,2)*Orb_eVec(:,2)',eye(2));
    Rel(ii,jj)=(H_eVec(:,aa)'*charup*H_eVec(:,aa))*2.37e-24*Vtp(ii,1)*Vt^2;
    %(H_eVec(:,aa)'*charup*H_eVec(:,aa))
    
    nuc_e=tensor([0 1],Orb_eVec(:,1)',[0 1]);
    nuc_g=tensor([1 0],Orb_eVec(:,1)',[0 1]);
    ele_e=tensor([0 1],Orb_eVec(:,1)',[1 0]);
    char_e=tensor([1 0],Orb_eVec(:,2)',[0 1]);
    nucfrac(ii,1)=abs(nuc_e*H_eVec(:,aa))^2;
    flipfrac(ii,1)=abs(nuc_e*H_eVec(:,cc))^2;
    %charfrac(ii,1)=abs(char_e*H_eVec(:,aa))^2;
    
    gB=ge*Bac/4;    
    deltaB=(Eud-Edd)-nuB(jj);
    deltaSO=Vtp(ii,1)-(gn*B0+nuB(jj)+Ap/2+deltaB);
    gSO=A/4*Vt/Vtp(ii,1);
    plus=(deltaSO+sqrt(deltaSO^2+(2*gSO)^2))/2/gSO;
    minus=(deltaSO-sqrt(deltaSO^2+(2*gSO)^2))/2/gSO;
    alpha(ii,1)=minus/sqrt(minus^2+1);
    beta(ii,1)=plus/sqrt(plus^2+1);
    DSO=deltaSO/2*(sqrt(1+(2*gSO/deltaSO)^2)-1);
    %(deltaB-DSO)/gB
    
    plusup=((deltaB+Ap)+sqrt((deltaB+Ap)^2+(2*gB)^2))/2/gB;
    minusup=((deltaB+Ap)-sqrt((deltaB+Ap)^2+(2*gB)^2))/2/gB;
    alphaup(ii,1)=minusup/sqrt(minusup^2+1);
    betaup(ii,1)=plusup/sqrt(plusup^2+1);
    
    plusdown1=((deltaB-DSO)+sqrt((deltaB-DSO)^2+(2*gB*beta(ii,1))^2))/2/gB/beta(ii,1);
    minusdown1=((deltaB-DSO)-sqrt((deltaB-DSO)^2+(2*gB*beta(ii,1))^2))/2/gB/beta(ii,1);
    alphadown1(ii,1)=minusdown1/sqrt(minusdown1^2+1);
    betadown1(ii,1)=plusdown1/sqrt(plusdown1^2+1);
    
    plusdown2=((deltaB+deltaSO+DSO)+sqrt((deltaB+deltaSO+DSO)^2+(-2*gB*alpha(ii,1))^2))/2/gB/(-alpha(ii,1));
    minusdown2=((deltaB+deltaSO+DSO)-sqrt((deltaB+deltaSO+DSO)^2+(-2*gB*alpha(ii,1))^2))/2/gB/(-alpha(ii,1));
    alphadown2(ii,1)=minusdown2/sqrt(minusdown2^2+1);
    betadown2(ii,1)=plusdown2/sqrt(plusdown2^2+1);
    %alphadown2.^2+betadown2.^2
    
    RelTheor(ii,jj)=(min(abs(betadown1(ii,1)),abs(alphadown1(ii,1)))*alpha(ii,1)+alphadown2(ii,1)*beta(ii,1))^2*2.37e-24*Vtp(ii,1)*Vt^2;
    
    hyper(ii,1) = Edd-Edu + nuB(jj);
    ensrot(ii,1) = gn*B0 + Ap/2 + (deltaB+Ap)/2*(sqrt(1+(2*gB/(deltaB+Ap))^2)-1) - (deltaB-DSO)/2*(sqrt(1+(2*gB*beta(ii,1)/(deltaB-DSO))^2)-1) - (deltaB+deltaSO+DSO)/2*(sqrt(1+(2*gB*(-alpha(ii,1))/(deltaB+deltaSO+DSO))^2)-1);
    %effrot(ii,1) = gn*B0 + sign(abs(betaup)-abs(alphaup))*sqrt((Ap+deltaB)^2+(2*gB)^2)/2 + sign(abs(betadown)-abs(alphadown))*sqrt(deltaB^2+(2*gB)^2)/2;
    
    gE=e*Eac*2*L/4/h*Vt/Vtp(ii,1);
    
    gEns(ii,jj)=-gE*betaup(ii,1)*(min(abs(betadown1(ii,1)),abs(alphadown1(ii,1)))*alpha(ii,1)+alphadown2(ii,1)*beta(ii,1));%
    
    nucfractheory(ii,1)=(max(abs(betadown1(ii,1)),abs(alphadown1(ii,1)))*betadown2(ii,1))^2;
    
    % find gate time
    nuE=ggB0(ii,1);
    
    H_expm=eye(8);
    for tt=1:T_coarse
        H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
        H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
    end
    
    tmax=(A/4*Vt/Vtp(ii,1))/(28e9*Bac/4)/(e*Eac*2*L/h/4*Vt/Vtp(ii,1))*100;
    tdim1=tmax*nuE;
    H_expm1=H_expm^(tdim1/tdim);
    
    Psit(:,1)=H_eVec(:,aa);
    for t=1:tdim
        Psit(:,t+1)=H_expm1*Psit(:,t);
        PsitB(t+1)=abs(Psit(:,t+1)'*nuc_g')^2;
    end
    
    AA=fft(PsitB);
    [pk,loc] = max(abs(AA(2:tdim/2)));
    nuRabi(ii,jj)=loc/tmax;
    
end
jj
toc;
end

%%
figure(23458)
subplot(2,2,1)
imagesc(Edc/1e3,nuB,nuRabi'/2/Eac);axis xy;title('Dipole numerics');colorbar
subplot(2,2,2)
imagesc(Edc/1e3,nuB,gEns'/Eac);axis xy;title('Dipole theory')
subplot(2,2,3)
imagesc(Edc/1e3,nuB,Rel');axis xy;title('Relaxation numerics');colorbar
subplot(2,2,4)
imagesc(Edc/1e3,nuB,RelTheor');axis xy;title('Relaxation theory')

%%
figure(23458)
subplot(2,2,1)
hold on
plot(Edc/1e3,nuB+deltaB, 'g', 'Linewidth', 1);%plot(Edc/1e3,nucfrac,Edc/1e3,nucfractheory)
plot(Edc/1e3,eo-hyper+nuB, 'b', 'Linewidth', 1);%plot(Edc/1e3,nucfrac,Edc/1e3,nucfractheory)
plot(Edc/1e3,edressed, 'w', 'Linewidth', 1);%plot(Edc/1e3,nucfrac,Edc/1e3,nucfractheory)
hold off

%% Noise nucleus rotating frame, resonant charge transition, sweep Edc - LEVEL DIAGRAM IN ROTATING FRAME - charge perturbs drive
clearvars
addpath Z:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

L=7.5e-9; %half interdot distance
A=117e6;
B0=1;%
Bac=0.6e-3;%0.477e-3;%
dg=-0.002; %relative difference in g-factor at interface

Eac = 32*0.001; % Electric drive
tdim=2000;
T_coarse=200;

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% Basis: nucleus x charge x electron
% 01) u,e,u
% 02) u,e,d
% 03) u,g,u
% 04) u,g,d
% 05) d,e,u
% 06) d,e,d
% 07) d,g,u
% 08) d,g,d

%Edc=570;

Vt=sqrt(((28e9*0.999+17.2e6)*B0)^2+A^2/4)+3e6;%
nuB=28e9*0.999*B0-A/4-2e6;
H_Znuc=tensor((-17.2e6*B0-nuB)/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z)+tensor(identity,identity,(-nuB)/2*sigma_z);
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_ESRnuc=tensor(-17.2e6*Bac/2/2*sigma_x,identity,identity);
H_ESRel=tensor(identity,identity-(identity/2+sigma_x/2)*0.002*0,28e9*Bac/2/2*sigma_x);
%[H_eVec,H_eVal]=eig(H_Znuc+H_tunel+H_Zel+H_Hyper+H_ESRnuc+H_ESRel);

% nuc_e=tensor([0 1],[0 1],[0 1]);
% nuc_g=tensor([1 0],[0 1],[0 1]);
% tic;
nuc_e=tensor([0 1],[0 1],[0 1]);
nuc_g=tensor([1 0],[0 1],[0 1]);
ele_e=tensor([0 1],[0 1],[1 0]);

Edcdim=200;
Edcmin=-2e4*B0/2;
Edcmax=2e4*B0/2;
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;

ggB0=zeros(Edcdim,1);
gnB0=zeros(Edcdim,1);
effrot=zeros(Edcdim,1);
ensrot=zeros(Edcdim,1);
hyper=zeros(Edcdim,1);
gEns=zeros(Edcdim,1);
Vtp=zeros(Edcdim,1);

nuRabi=zeros(Edcdim,1);
Psit=zeros(8,tdim+1);
PsitB=zeros(1,tdim+1);

for ii=1:Edcdim
    tic;
    H_Edc=-tensor(identity,e*Edc(ii)*L/h*sigma_x,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
    [H_eVec,H_eVal]=eig(H_0);
    
    fid1=0;
    fid2=0;
    fid3=0;
    aa=0;
    bb=0;
    cc=0;
    for i=1:8
        AA=abs(nuc_e*H_eVec(:,i))^2;
        BB=abs(nuc_g*H_eVec(:,i))^2;
        CC=abs(ele_e*H_eVec(:,i))^2;
        if AA>fid1
            fid1=AA;
            aa=i;
        end
        if BB>fid2
            fid2=BB;
            bb=i;
        end
        if CC>fid3
            fid3=CC;
            cc=i;
        end
    end
    
    ggB0(ii,1)=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
    gnB0(ii,1)=abs(H_eVal(cc,cc)-H_eVal(bb,bb));
    
    Vtp(ii,1) = sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2);
    Ap = A/2*(1-((e*2*L*Edc(ii)/h)/Vtp(ii,1)));
    ge = 28e9*(1+(1+((e*2*L*Edc(ii)/h)/Vtp(ii,1)))*dg/2);
    g_ffo=A/4*Vt/Vtp(ii,1);
    gn = 17.2e6;
    %electron,nucleus
    Euu =  ((ge - gn)*B0 + Ap/2)/2;
    Edd =  (-1*(ge - gn)*B0 + Ap/2)/2;
    Edu =  (-1*sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    Eud =  (sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    
    gB=ge*Bac/4;
    deltaB=(Eud-Edd)-nuB;
    plusdown=(deltaB+sqrt(deltaB^2+(2*gB)^2))/2/gB;
    minusdown=(deltaB-sqrt(deltaB^2+(2*gB)^2))/2/gB;
    plusup=((deltaB+Ap)+sqrt((deltaB+Ap)^2+(2*gB)^2))/2/gB;
    minusup=((deltaB+Ap)-sqrt((deltaB+Ap)^2+(2*gB)^2))/2/gB;
    alphadown=minusdown/sqrt(minusdown^2+1);
    betadown=plusdown/sqrt(plusdown^2+1);
    alphaup=minusup/sqrt(minusup^2+1);
    betaup=plusup/sqrt(plusup^2+1);
    %alphaup^2+betaup^2
    %alphadown^2+betadown^2
    
    hyper(ii,1) = Edd-Edu + nuB;
    ensrot(ii,1) = gn*B0 + sign(abs(betaup)-abs(alphaup))*sqrt((Ap+deltaB)^2+(2*gB)^2)/2 + sign(abs(alphadown)-abs(betadown))*sqrt(deltaB^2+(2*gB)^2)/2;
    deltaSOns=Vtp(ii,1)-ensrot(ii,1)-nuB;
    gsons=A/4*Vt/Vtp(ii,1)*max(abs(betaup),abs(alphaup))*min(abs(betadown),abs(alphadown));
    ensrot(ii,1) = ensrot(ii,1) - deltaSOns/2*(sqrt(1+(2*gsons/deltaSOns)^2)-1);
    effrot(ii,1) = gn*B0 + sign(abs(betaup)-abs(alphaup))*sqrt((Ap+deltaB)^2+(2*gB)^2)/2 + sign(abs(betadown)-abs(alphadown))*sqrt(deltaB^2+(2*gB)^2)/2;
    deltaSOff=Vtp(ii,1)-effrot(ii,1)-nuB;
    gsoff=A/4*Vt/Vtp(ii,1)*max(abs(betaup),abs(alphaup))*max(abs(betadown),abs(alphadown));
    effrot(ii,1) = effrot(ii,1) - deltaSOff/2*(sqrt(1+(2*gsoff/deltaSOff)^2)-1);
    
    gE=e*Eac*2*L/4/h*Vt/Vtp(ii,1);
    deltaE=Vtp(ii,1)-ggB0(ii,1);
    
    gEns(ii,1)=gE/deltaSOns*gsons;%gE/deltaE*gsons;
    
    % find gate time
    nuE=ggB0(ii,1);
    
    H_expm=eye(8);
    for tt=1:T_coarse
        H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
        H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
    end
    
    tmax=(A/4*Vt/Vtp(ii,1))/(28e9*Bac/4)/(e*Eac*2*L/h/4*Vt/Vtp(ii,1))*100;
    tdim1=tmax*nuE;
    H_expm1=H_expm^(tdim1/tdim);
    
    Psit(:,1)=H_eVec(:,aa);
    for t=1:tdim
        Psit(:,t+1)=H_expm1*Psit(:,t);
        PsitB(t+1)=abs(Psit(:,t+1)'*nuc_g')^2;
    end
    
    AA=fft(PsitB);
    [pk,loc] = max(abs(AA(2:tdim/2)));
    nuRabi(ii,1)=loc/tmax;
    
    ii
    toc;
end
%%
figure(111)
subplot(1,2,1)
plot(Edc,hyper-nuB,Edc,ggB0-nuB,Edc,gnB0-nuB,Edc,effrot,Edc,ensrot,Edc,Vtp-nuB);%
legend('bare','numerics','numerics','theory','theory','charge')%
ylim([min(min(ggB0),min(gnB0))-nuB, max(max(ggB0),max(gnB0))-nuB])
subplot(1,2,2)
plot(Edc,nuRabi/2/Eac,Edc,gEns/Eac)%
%ylim([0,2e5]);

%% Noise nucleus rotating frame, resonant charge transition, sweep Edc
clearvars
addpath Z:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

L=7.5e-9; %half interdot distance
A=117e6;
B0=1.2;%
Bac=0.5e-3;%0.477e-3;%
dg=-0.002; %relative difference in g-factor at interface

Eac = 32*0.001; % Electric drive
tdim=2000;
T_coarse=200;

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% Basis: nucleus x charge x electron
% 01) u,e,u
% 02) u,e,d
% 03) u,g,u
% 04) u,g,d
% 05) d,e,u
% 06) d,e,d
% 07) d,g,u
% 08) d,g,d

%Edc=570;

Vt=sqrt(((28e9*0.999+17.2e6)*B0)^2+A^2/4)+3e6;%
nuB=28e9*0.999*B0-A/4-2.5e6;
H_Znuc=tensor((-17.2e6*B0-nuB)/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z)+tensor(identity,identity,(-nuB)/2*sigma_z);
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_ESRnuc=tensor(-17.2e6*Bac/2/2*sigma_x,identity,identity);
H_ESRel=tensor(identity,identity-(identity/2+sigma_x/2)*0.002*0,28e9*Bac/2/2*sigma_x);
%[H_eVec,H_eVal]=eig(H_Znuc+H_tunel+H_Zel+H_Hyper+H_ESRnuc+H_ESRel);

% nuc_e=tensor([0 1],[0 1],[0 1]);
% nuc_g=tensor([1 0],[0 1],[0 1]);
% tic;
nuc_e=tensor([0 1],[0 1],[0 1]);
nuc_g=tensor([1 0],[0 1],[0 1]);
ele_e=tensor([0 1],[0 1],[1 0]);

Edcdim=100;
Edcmin=-2e4*B0/10;
Edcmax=2e4*B0/10;
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;

ggB0=zeros(Edcdim,1);
gnB0=zeros(Edcdim,1);
stark_toe1=zeros(Edcdim,1);
stark_toe2=zeros(Edcdim,1);
hyper=zeros(Edcdim,1);
gEns=zeros(Edcdim,1);

nuRabi=zeros(Edcdim,1);
Psit=zeros(8,tdim+1);
PsitB=zeros(1,tdim+1);

for ii=1:Edcdim
    tic;
    H_Edc=-tensor(identity,e*Edc(ii)*L/h*sigma_x,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
    [H_eVec,H_eVal]=eig(H_0);
    
    fid1=0;
    fid2=0;
    fid3=0;
    aa=0;
    bb=0;
    cc=0;
    for i=1:8
        AA=abs(nuc_e*H_eVec(:,i))^2;
        BB=abs(nuc_g*H_eVec(:,i))^2;
        CC=abs(ele_e*H_eVec(:,i))^2;
        if AA>fid1
            fid1=AA;
            aa=i;
        end
        if BB>fid2
            fid2=BB;
            bb=i;
        end
        if CC>fid3
            fid3=CC;
            cc=i;
        end
    end
    
    ggB0(ii,1)=abs(H_eVal(aa,aa)-H_eVal(bb,bb));
    gnB0(ii,1)=abs(H_eVal(cc,cc)-H_eVal(bb,bb));
    
    Vtp = sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2);
    Ap = A/2*(1-((e*2*L*Edc(ii)/h)/Vtp));
    ge = 28e9*(1+(1+((e*2*L*Edc(ii)/h)/Vtp))*dg/2);
    g_ffo=A/4*Vt/Vtp;
    gn = 17.2e6;
    %electron,nucleus
    Euu =  ((ge - gn)*B0 + Ap/2)/2;
    Edd =  (-1*(ge - gn)*B0 + Ap/2)/2;
    Edu =  (-1*sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
    Eud =  (sqrt(((ge + gn)*B0)^2 + (Ap)^2) - Ap/2)/2;
 
    Delta3=(Euu-Edu)-nuB;%nuB-(ge*B0+Ap/2);
    gnu=ge*Bac/4;
    %Vtp=Vtp+gnu^2/Delta3;
    
    SO_level1=((Eud-Edu)+Vtp)/2-sqrt(((Eud-Edu)-Vtp)^2+(2*g_ffo)^2)/2;%(sqrt(((ge+17.2e6)*B0)^2+Ap^2)+Vtp)/2-sqrt((sqrt(((ge+17.2e6)*B0)^2+Ap^2)-Vtp)^2+(2*g_ffo)^2)/2;
    SO_level2=((Eud-Edu)+Vtp)/2+sqrt(((Eud-Edu)-Vtp)^2+(2*g_ffo)^2)/2;%(sqrt(((ge+17.2e6)*B0)^2+Ap^2)+Vtp)/2+sqrt((sqrt(((ge+17.2e6)*B0)^2+Ap^2)-Vtp)^2+(2*g_ffo)^2)/2;
    Delta1=nuB+(Edd-Edu)-SO_level1;%nuB-(SO_level1-(Ap/2+17.2e6*B0));
    Delta2=nuB+(Edd-Edu)-SO_level2;%nuB-(SO_level2-(Ap/2+17.2e6*B0));
    %Sigma1=nuB+(Edd-Edu)+SO_level1;%nuB-(SO_level1-(Ap/2+17.2e6*B0));
    %Sigma2=nuB+(Edd-Edu)+SO_level2;%nuB-(SO_level2-(Ap/2+17.2e6*B0));

    aaa=(((Eud-Edu)-Vtp)-sqrt(((Eud-Edu)-Vtp)^2+(2*g_ffo)^2))/2/g_ffo;
    bbb=(((Eud-Edu)-Vtp)+sqrt(((Eud-Edu)-Vtp)^2+(2*g_ffo)^2))/2/g_ffo;
    alpha=aaa/sqrt(aaa^2+1);
    beta=bbb/sqrt(bbb^2+1);
    gnu1=gnu*alpha;
    gnu2=gnu*beta;
    hyper(ii,1) = Edd-Edu + nuB;
    stark_toe2(ii,1)  = hyper(ii,1) + Delta3*(sqrt(1+(2*gnu)^2/Delta3^2)-1)/2 + Delta1*(sqrt(1+(2*gnu1)^2/Delta1^2)-1)/2 + Delta2*(sqrt(1+(2*gnu2)^2/Delta2^2)-1)/2;% - (ge*Bac/4)^2/(nuB-(ge*B0+Ap/2))
    stark_toe1(ii,1) = SO_level1 + Delta3*(sqrt(1+(2*gnu)^2/Delta3^2)-1)/2 - Delta1*(sqrt(1+(2*gnu1)^2/Delta1^2)-1)/2 - 0*Delta2*(sqrt(1+(2*gnu2)^2/Delta2^2)-1)/2;% + Delta1*(sqrt(1+(2*gnu1)^2/Delta1^2)-1)/2%- Delta3*(sqrt(1+(2*gnu)^2/Delta3^2)-1)/2
%     Delta3*(sqrt(1+(2*gnu)^2/Delta3^2)-1)/2
%     Delta1*(sqrt(1+(2*gnu1)^2/Delta1^2)-1)/2
%     Delta2*(sqrt(1+(2*gnu2)^2/Delta2^2)-1)/2
    
    gE=e*Eac*2*L/4/h*Vt/Vtp;
    gEns(ii,1)=-gE*gnu*alpha*beta*abs(1/Delta2);%-1/Sigma2+1/Sigma1
    gE*beta/Delta1
    gnu*alpha/Delta1
    beta
    %aaa/sqrt(aaa^2+1)*bbb/sqrt(bbb^2+1)
    
    % find gate time
    nuE=ggB0(ii,1);
    
    H_expm=eye(8);
    for tt=1:T_coarse
        H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
        H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
    end
    
    tmax=(A/4*Vt/Vtp)/(28e9*Bac/4)/(e*Eac*2*L/h/4*Vt/Vtp)*100;
    tdim1=tmax*nuE;
    H_expm1=H_expm^(tdim1/tdim);
    
    Psit(:,1)=H_eVec(:,aa);
    for t=1:tdim
        Psit(:,t+1)=H_expm1*Psit(:,t);
        PsitB(t+1)=abs(Psit(:,t+1)'*nuc_g')^2;
    end
    
    AA=fft(PsitB);
    [pk,loc] = max(abs(AA(2:tdim/2)));
    nuRabi(ii,1)=loc/tmax;
    
    ii
    toc;
end
%%
figure(111)
subplot(1,2,1)
plot(Edc,hyper-nuB,Edc,ggB0-nuB,Edc,stark_toe2-nuB);%,Edc,gnB0-nuB,Edc,stark_toe1-nuB
legend('bare','numerics','theory')%,'numerics','theory'
subplot(1,2,2)
plot(Edc,nuRabi/2/Eac,Edc,gEns/Eac)%
ylim([0,2e5]);

%% Electric drive
Edc=177;

H_Edc=-tensor(identity,e*Edc*L/h*sigma_x,identity);
H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
[H_eVec,H_eVal]=eig(H_0);

fid1=0;
fid2=0;
aa=0;
bb=0;
for i=1:8
    AA=abs(nuc_e*H_eVec(:,i))^2;
    BB=abs(nuc_g*H_eVec(:,i))^2;
    if AA>fid1
        fid1=AA;
        aa=i;
    end
    if BB>fid2
        fid2=BB;
        bb=i;
    end
end
nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb));%-2.2e5

tdim=1000;
T_coarse=200;
H_expm=eye(8);
for tt=1:T_coarse
    H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
    H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
end

eo=sqrt(Vt^2+(e*Edc*2*L/h)^2);
tmax=(A/4*Vt/eo)/(28e9*Bac/4)/(e*Eac*2*L/h/4*Vt/eo)*40;
%tmax=round(tmax*nuE/tdim)*tdim/nuE;%makes sure time evolution finishes at a nuE 2pi period
tdim1=tmax*nuE;
H_expm1=H_expm^(tdim1/tdim);

Psit=zeros(8,tdim+1);
Psit_diag=zeros(8,tdim+1);
PsitB=zeros(1,tdim+1);

Psit(:,1)=H_eVec(:,bb);
Psit_diag(:,1)=ctranspose(H_eVec)*Psit(:,1);
PsitB(:,1)=abs(Psit(:,1)'*nuc_g')^2;

for t=1:tdim
    Psit(:,t+1)=H_expm1*Psit(:,t);
    Psit_diag(:,t+1)=ctranspose(H_eVec)*Psit(:,t);
    PsitB(t+1)=abs(Psit(:,t+1)'*nuc_g')^2;
end

figure(222);%
subplot(7,1,1);
plot(0:tmax/tdim:tmax,abs(Psit(:,:).*Psit(:,:)).');%legend('1','2','3','4','5','6','7','8')
subplot(7,1,2);
plot(0:tmax/tdim:tmax,abs(Psit_diag(:,:).*Psit_diag(:,:)).');
subplot(7,1,3);
imagesc(abs(Psit(:,:).*Psit(:,:)).');%_diag
subplot(7,1,4);
%plot(0:tmax/tdim:tmax,PsitB);
AA=fft(PsitB);
plot(abs(AA(2:tdim/2)));
[pk,loc] = max(abs(AA(2:tdim/2)));
nuRabi=loc/tmax;
nuRabi/Eac/2
tmaxpi2=1/nuRabi;

tmax=tmaxpi2;
tdim2=fix(tmax*nuE);
H_expm2=H_expm^(tdim2/tdim);

for t=1:tdim
    Psit(:,t+1)=H_expm2*Psit(:,t);
    Psit_diag(:,t+1)=ctranspose(H_eVec)*Psit(:,t);
end

tmaxpi2=fix(tmaxpi2*nuE)/nuE;
tdim2=tmaxpi2*nuE;

H_expm_extra=eye(8);
for t=1:round((1/4/nuRabi-tmaxpi2)*nuE*T_coarse)
    H=H_0+cos(2*pi*t/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
	H_expm_extra=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm_extra;
end

Psipi2=H_expm_extra*H_expm^tdim2*Psit(:,1);
%fidpi2=max([abs(dot(H_eVec(:,1)-exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2))^2/2 abs(dot(H_eVec(:,1)+exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2))^2/2])

figure(222);
subplot(7,1,5);
plot(0:tmax/tdim:tmax,abs(Psit(:,:).*Psit(:,:)).');
subplot(7,1,6);
plot(0:tmax/tdim:tmax,abs(Psit_diag(:,:).*Psit_diag(:,:)).');
subplot(7,1,7);
imagesc(abs(Psit(:,:).*Psit(:,:)).');%_diag
%fidpi2=max([abs(dot(H_eVec(:,1)+i*H_eVec(:,4),Psipi2))^2/2 abs(dot(H_eVec(:,1)-i*H_eVec(:,4),Psipi2))^2/2])%abs(dot(tensor([0 1],[0 1],[1 0])+i*tensor([1 0],[0 1],[0 1]),Psit(:,tdim+1)))^2/2

%% Electric drive - bare levels?
Edc=1000;

H_Edc=-tensor(identity,e*Edc*L/h*sigma_x,identity);

H_Znuc=tensor((-17.2e6*B0)/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z);
H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;%+H_ESRnuc+H_ESRel;
[H_eVec,H_eVal]=eig(H_0);

fid1=0;
fid2=0;
aa=0;
bb=0;
for i=1:8
    AA=abs(nuc_e*H_eVec(:,i))^2;
    BB=abs(nuc_g*H_eVec(:,i))^2;
    if AA>fid1
        fid1=AA;
        aa=i;
    end
    if BB>fid2
        fid2=BB;
        bb=i;
    end
end
nuE=abs(H_eVal(aa,aa)-H_eVal(bb,bb))+nuB

H_Znuc=tensor((-17.2e6*B0-nuB)/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z)+tensor(identity,identity,(-nuB)/2*sigma_z);
H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;

tdim=1000;
T_coarse=200;
H_expm=eye(8);
for tt=1:T_coarse
    H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
    H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
end

eo=sqrt(Vt^2+(e*Edc*2*L/h)^2);
tmax=(A/4*Vt/eo)/(28e9*Bac/4)/(e*Eac*2*L/h/4*Vt/eo)*40;
%tmax=round(tmax*nuE/tdim)*tdim/nuE;%makes sure time evolution finishes at a nuE 2pi period
tdim1=tmax*nuE;
H_expm1=H_expm^(tdim1/tdim);

Psit=zeros(8,tdim+1);
Psit_diag=zeros(8,tdim+1);
PsitB=zeros(1,tdim+1);

Psit(:,1)=H_eVec(:,aa);
Psit_diag(:,1)=ctranspose(H_eVec)*Psit(:,1);
PsitB(:,1)=abs(Psit(:,1)'*nuc_g')^2;

for t=1:tdim
    Psit(:,t+1)=H_expm1*Psit(:,t);
    Psit_diag(:,t+1)=ctranspose(H_eVec)*Psit(:,t);
    PsitB(t+1)=abs(Psit(:,t+1)'*nuc_g')^2;
end

figure(222);%
subplot(7,1,1);
plot(0:tmax/tdim:tmax,abs(Psit(:,:).*Psit(:,:)).');%legend('1','2','3','4','5','6','7','8')
subplot(7,1,2);
plot(0:tmax/tdim:tmax,abs(Psit_diag(:,:).*Psit_diag(:,:)).');
subplot(7,1,3);
imagesc(abs(Psit(:,:).*Psit(:,:)).');%_diag
subplot(7,1,4);
%plot(0:tmax/tdim:tmax,PsitB);
AA=fft(PsitB);
plot(abs(AA(2:tdim/2)));
[pk,loc] = max(abs(AA(2:tdim/2)));
nuRabi=loc/tmax
tmaxpi2=1/nuRabi;

tmax=tmaxpi2;
tdim2=fix(tmax*nuE);
H_expm2=H_expm^(tdim2/tdim);

for t=1:tdim
    Psit(:,t+1)=H_expm2*Psit(:,t);
    Psit_diag(:,t+1)=ctranspose(H_eVec)*Psit(:,t);
end

tmaxpi2=fix(tmaxpi2*nuE)/nuE;
tdim2=tmaxpi2*nuE;

H_expm_extra=eye(8);
for t=1:round((1/4/nuRabi-tmaxpi2)*nuE*T_coarse)
    H=H_0+cos(2*pi*t/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
	H_expm_extra=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm_extra;
end

Psipi2=H_expm_extra*H_expm^tdim2*Psit(:,1);
%fidpi2=max([abs(dot(H_eVec(:,1)-exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2))^2/2 abs(dot(H_eVec(:,1)+exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2))^2/2])

figure(222);
subplot(7,1,5);
plot(0:tmax/tdim:tmax,abs(Psit(:,:).*Psit(:,:)).');
subplot(7,1,6);
plot(0:tmax/tdim:tmax,abs(Psit_diag(:,:).*Psit_diag(:,:)).');
subplot(7,1,7);
imagesc(abs(Psit(:,:).*Psit(:,:)).');%_diag
%fidpi2=max([abs(dot(H_eVec(:,1)+i*H_eVec(:,4),Psipi2))^2/2 abs(dot(H_eVec(:,1)-i*H_eVec(:,4),Psipi2))^2/2])%abs(dot(tensor([0 1],[0 1],[1 0])+i*tensor([1 0],[0 1],[0 1]),Psit(:,tdim+1)))^2/2


%% sweep Edc noise and average
Edcnoisedim=3;%it has to be an odd number to include the zero noise point to the average that is supposed in the code!
Edcnoise=-1:2/(Edcnoisedim-1):1;
RMSnoise=100;
Edcnoisemax=RMSnoise/sqrt(sum(Edcnoise.^2)/Edcnoisedim)
Edcnoise=(-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax;

Dephasing=0;
fidpi2noise=0;
fidpi2noiseExtra=0;

for kk=1:Edcnoisedim
    H_Edc=tensor(identity,e*Edcnoise(kk)*L/h*sigma_x,identity);
    H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
    [H_eVec,H_eVal]=eig(H_0);
    nuEnoise=abs(H_eVal(1,1)-H_eVal(4,4));
    Psit=H_eVec(:,1);
    
    H_expm=eye(8);
    for tt=1:T_coarse
        H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
        H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
    end
    H_expm_extra=eye(8);
    for t=1:round((pi/2/CC(2)-tmaxpi2)*nuE*T_coarse)
        H=H_0+cos(2*pi*t/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
        H_expm_extra=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm_extra;
    end
    
    Dephasing=Dephasing+abs(nuE-nuEnoise)/Edcnoisedim;
    
    Psipi2noise=H_expm_extra*H_expm^tdim2*H_eVec(:,1);%Psipi2noise=H_expm^tdim1*Psit(:,1);
    fidpi2noise=fidpi2noise+max([abs(dot(H_eVec(:,1)-exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2noise))^2/2 abs(dot(H_eVec(:,1)+exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2noise))^2/2])/Edcnoisedim;
    fidpi2noiseExtra=fidpi2noiseExtra+abs(dot(Psipi2,Psipi2noise))^2/Edcnoisedim;
end
Dephasing
fidpi2noise
%% 1-qubit gate fidelity assuming eigenstates (nuclear-spin qubit)
clear all;
addpath Q:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

L=7.5e-9; %half interdot distance
A=117e6;
B0=1;%1
dg=-0.002; %relative difference in g-factor at interface
Vt=sqrt(((28e9*0.999+17.2e6)*B0)^2+A^2/4)+5e6;
nuB=28e9*B0*0.999-A/4-0e6;%9.82e+06;%28e9*B0-A/4;%-1.97e+07

Eac = 32*0.025;

Edcnoisedim=21;
Edcnoisemax=170;

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);
 
H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Znuc=tensor((-17.2e6*B0-nuB)/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,(28e9*B0)/2*sigma_z)+tensor(identity,identity,(-nuB)/2*sigma_z);
%H_0=H_Znuc+H_tunel+H_Zel+H_Hyper;%+H_Edc+H_ESRnuc+H_ESRel;

tdim=2000;
%Psit=zeros(8,tdim);
%Psit_diag=zeros(8,tdim);

%  Vt & Edc sweep
dimBac=15*4*10/20*6;
%ii=1:dimBac;
Bacmin=0.49e-3;
Bacmax=0.73e-3;
Bac=Bacmin:(Bacmax-Bacmin)/(dimBac-1):Bacmax;
dimEdc=14*4*5/10*4;
%jj=1:dimEdc;
Edcmin=-930;
Edcmax=-230;
Edc=Edcmin:(Edcmax-Edcmin)/(dimEdc-1):Edcmax;
fidpi2=zeros(dimBac,dimEdc);
fidpi2noise=zeros(dimBac,dimEdc);
fidpi2noiseExtra=zeros(dimBac,dimEdc);
pi2time=zeros(dimBac,dimEdc);
Dephasing=zeros(dimBac,dimEdc);
RabiJitter=zeros(dimBac,dimEdc);

T_coarse=300;
%Psit(:,1)=tensor([1 0],[0 1],[0 1]);%H_eVec(:,1);%tensor([0 1],[0 1],[1 0]);%%-i*H_eVec(:,3))/sqrt(2)
%Psit_diag(:,1)=ctranspose(H_eVec)*Psit(:,1);

%delete(gcp);parpool(28)

tic;
for ii=1:dimBac
    ii
    H_ESRnuc=tensor(-17.2e6*Bac(ii)/2/2*sigma_x,identity,identity);
    H_ESRel=tensor(identity,identity-(identity/2+sigma_x/2)*0.002*0,28e9*Bac(ii)/2/2*sigma_x);
    parfor jj=1:dimEdc %parfor
        H_Edc=tensor(identity,e*Edc(jj)*L/h*sigma_x,identity);
        
        H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
        
        [H_eVec,H_eVal]=eig(H_0);
        
        nuE=abs(H_eVal(1,1)-H_eVal(4,4));%nuE=(28e9+17.2e6)*B0*0.9985;%+5e6;
        
        H_expm=eye(8);
        for tt=1:T_coarse
            H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
            H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
        end
        
        eo=sqrt(Vt^2+(e*Edc(jj)*2*L/h)^2);
        tmax=(A/4*Vt/eo)/(28e9*Bac(ii)/4)/(e*Eac*2*L/h/4*Vt/eo)*40;%500e-6;
        tmax=round(tmax*nuE/tdim)*tdim/nuE;%makes sure time evolution finishes at a nuE 2pi period
        tdim1=tmax*nuE;
        H_expm_dyn=H_expm^(tdim1/tdim);
        
        Psit=H_eVec(:,1);
        %Psit_diag(:,1)=ctranspose(H_eVec)*Psit(:,1);
        for t=1:tdim
            Psit(:,t+1)=H_expm_dyn*Psit(:,t);
            %Psit_diag(:,t+1)=ctranspose(H_eVec)*Psit(:,t);
        end

        % find Rabi frequency
%         AA=fft(abs(Psit(4,:)).^2);
%         [pks,locs] = findpeaks(abs(AA),'MinPeakProminence',30);
%         fRabi=locs(1);
%         if size(locs)==[3 1]
%             fRabi=min([abs(locs(1)-locs(2)), abs(locs(1)-locs(3)), abs(locs(3)-locs(2))]);
%         else
%             fRabi=min([abs(locs(1)-locs(2)), abs(locs(1)-0.5), abs(locs(2)-0.5), abs(1000.5-locs(1)), abs(1000.5-locs(2))]);
%         end
        
        % find state after pi/2 pulse
%         tmax=tmax/4/fRabi;
        shft=min(abs(Psit(4,:)).^2);
        CC=coeffvalues(fit((0:tmax/tdim:tmax).',(abs(Psit(4,:)).^2).'-(1+shft)/2,'sin1'));
        %fRabi=CC(2)
        tmaxpi2=pi/2/CC(2);%tmax/4/fRabi;%
        pi2time(ii,jj)=tmaxpi2;
        tmaxpi2=fix(tmaxpi2*nuE)/nuE;
        tdim2=tmaxpi2*nuE;
                
        H_expm_extra=eye(8);
        for t=1:round((pi/2/CC(2)-tmaxpi2)*nuE*T_coarse)
            H=H_0+cos(2*pi*t/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
            H_expm_extra=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm_extra;
        end
        
        Psipi2=H_expm_extra*H_expm^tdim2*Psit(:,1);
        fidpi2(ii,jj)=max([abs(dot(H_eVec(:,1)-exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2))^2/2 abs(dot(H_eVec(:,1)+exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2))^2/2]);
        
        % sweep Edc noise and average
        Edcnoise=Edc(jj)+((-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax);
        
        for kk=1:Edcnoisedim
            H_Edc=tensor(identity,e*Edcnoise(kk)*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc+H_ESRnuc+H_ESRel;
            [H_eVec,H_eVal]=eig(H_0);
            nuEnoise=abs(H_eVal(1,1)-H_eVal(4,4));
            
            H_expm=eye(8);
            for tt=1:T_coarse
                H=H_0+cos(2*pi*tt/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
                H_expm=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm;
            end
            H_expm_extra=eye(8);
            for t=1:round((pi/2/CC(2)-tmaxpi2)*nuE*T_coarse)
                H=H_0+cos(2*pi*t/T_coarse)*tensor(identity,e*Eac*L/h*sigma_x,identity);
                H_expm_extra=expm(-1i*2*pi*H/(nuE*T_coarse))*H_expm_extra;
            end
            
            % find new Rabi frequency
            Psit=H_eVec(:,1);
            H_expm_dyn=H_expm^(tdim1/tdim);
            for t=1:tdim
                Psit(:,t+1)=H_expm_dyn*Psit(:,t);
                %Psit_diag(:,t+1)=ctranspose(H_eVec)*Psit(:,t);
            end
            AA=fft(abs(Psit(4,:)).^2);
            [pks,locs] = findpeaks(abs(AA),'MinPeakProminence',10);
            if size(locs)==[3 1]
                fRabi=min([abs(locs(1)-locs(2)), abs(locs(1)-locs(3)), abs(locs(3)-locs(2))])
            else
                fRabi=min([abs(locs(1)-locs(2)), abs(locs(1)-0.5), abs(locs(2)-0.5), abs(1000.5-locs(1)), abs(1000.5-locs(2))])
            end
            tmaxpi2noise=tmax/4/fRabi;
            RabiJitter(ii,jj)=RabiJitter(ii,jj)+abs(tmaxpi2-tmaxpi2noise)/Edcnoisedim;

            Dephasing(ii,jj)=Dephasing(ii,jj)+abs(nuE-nuEnoise)/Edcnoisedim;
            
            Psipi2noise=H_expm_extra*H_expm^tdim2*H_eVec(:,1);%Psipi2noise=H_expm^tdim1*Psit(:,1);
            fidpi2noise(ii,jj)=fidpi2noise(ii,jj)+max([abs(dot(H_eVec(:,1)-exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2noise))^2/2 abs(dot(H_eVec(:,1)+exp(-1i*(2*pi*nuE*pi/2/CC(2)+pi/2))*H_eVec(:,4),Psipi2noise))^2/2])/Edcnoisedim;
            fidpi2noiseExtra(ii,jj)=fidpi2noiseExtra(ii,jj)+abs(dot(Psipi2,Psipi2noise))^2/Edcnoisedim;
        end
    end
end
toc;
    
%figure(gcf);
figure(333);
subplot(2,3,1);
imagesc(-Edc,Bac*1e3,log10(1-fidpi2));colorbar;%caxis([-3 -1])
xlabel('Edc(V/m)');ylabel('Bac(mT)');title('Intrinsic error');
subplot(2,3,2);
imagesc(-Edc,Bac*1e3,log10(1-fidpi2noiseExtra));colorbar;%caxis([-3 -1])
xlabel('Edc(V/m)');ylabel('Bac(mT)');title('Error added by noise');
subplot(2,3,3);
imagesc(-Edc,Bac*1e3,log10(1-fidpi2noise));colorbar;%caxis([-3 -1])
xlabel('Edc(V/m)');ylabel('Bac(mT)');title('Total error');
subplot(2,3,4);
imagesc(-Edc,Bac*1e3,pi2time*1e9);colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Bac(mT)');title('Gate time (ns)');
subplot(2,3,5);
imagesc(-Edc,Bac*1e3,RabiJitter./pi2time);colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Bac(mT)');title('Rabi jitter (ns)');
subplot(2,3,6);
imagesc(-Edc,Bac*1e3,log10(Dephasing));colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Bac(mT)');title('log[Dephasing rate (Hz)]');

