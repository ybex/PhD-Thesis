%% Dispersion sweep Edc 
% Figure 2a
clearvars
%addpath Z:\spin-QED\Matlab_functions

clear Edc Vtt ggB0
e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

L=7.5e-9; %half interdot distance
A=117e6; %contact hyperfine
B0=0.2*2;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+243e6; %double-dot tunel rate
dg=-0.002;%relative change in g-factor (delta_g/g)

Edcdim=1200;
Edcmin=-20000*B0;%332;%
Edcmax=20000*B0;%342;%
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;

Vtt=zeros(Edcdim,1);
fidVtt=zeros(Edcdim,1);
ggB0=zeros(Edcdim,1);
fidggB0=zeros(Edcdim,1);
dgEdc=zeros(Edcdim,1);
dStark=zeros(Edcdim,1);
dStarkVt=zeros(Edcdim,1);
dA=zeros(Edcdim,1);

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_tunel=tensor(identity,Vt/2*sigma_z,identity);
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

orb_e=tensor([1 0],[1 0],[0 1]);
flip_e=tensor([0 1],[0 1],[1 0]);
orbflip_g=tensor([1 0],[0 1],[0 1]);

tic;

for ii=1:Edcdim
    H_Edc=-tensor(identity,e*Edc(ii)*L/h*sigma_x,identity);
    H=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    
    fid1=0;
    fid2=0;
    fid3=0;
    aa=0;
    bb=0;
    cc=0;
    for i=1:8
        AA=abs(orb_e*H_eVec(:,i))^2;
        BB=abs(flip_e*H_eVec(:,i))^2;
        CC=abs(orbflip_g*H_eVec(:,i))^2;
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
    
    Vtt(ii,1)=abs(H_eVal(aa,aa)-H_eVal(cc,cc));
    fidVtt(ii,1)=fid1*fid3;
    ggB0(ii,1)=abs(H_eVal(bb,bb)-H_eVal(cc,cc));
    fidggB0(ii,1)=fid2*fid3;
    
    %fahd
    Vtp = sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2);
    gp = 1+(1/2+((e*2*L*Edc(ii)/h)/Vtp)/2)*dg;
    Ap = A/2*(1-((e*2*L*Edc(ii)/h)/Vtp));
    
    dA(ii,1) = sqrt(((28e9+17.2e6)*B0)^2+Ap^2);
    dgEdc(ii,1) = sqrt(((gp*28e9+17.2e6)*B0)^2+Ap^2);
    dStark(ii,1) = -(A/4*(Vt/Vtp))^2/(Vtp - (gp*28e9 + 17.2e6)*B0);
    dStarkVt(ii,1) = -(A/4)^2/(Vtp - (gp*28e9 + 17.2e6)*B0);
    
    
%     %
%     dgEdc(ii,1)=-e*2*L*Edc(ii)/h/sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2)/2*dg*28e9*B0;
%     dStark(ii,1)=-(A/4*Vt/sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2))^2/(sqrt(Vt^2+(e*2*L*Edc(ii)/h)^2)-(28e9+17.2e6)*B0);
end
toc;

figure(7771);
plot(Edc/1000,dgEdc/1e9,'Color',[.5 .5 .5],'LineWidth',1.5);
hold on
plot(Edc(1:1:1200)/1000,ggB0(1:1:1200)/1e9,'Color',[.5 .5 .5],'LineWidth',5);%,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerEdgeColor',[0 0 0],'LineWidth',1);
plot(Edc/1000,dgEdc/1e9+dStark/1e9,'Color',[0 0 0],'LineWidth',1.5);
hold off
legend('A(E),g(E)','Hamiltonian','A(E),g(E),disp(E)','Location','southwest');
%ylim([27.962e9 28.017e9]/1e9);
xlim([Edcmin Edcmax]/1000);
%set(gca, 'Position', [0.12 0.1 0.8 0.54],'YTick',[27.97 27.98 27.99 28.00 28.01])

%% Figure 2a plotting
FigHandle = figure(777);
set(FigHandle, 'Position', [600, 100, 450, 900]);

subplot(2,1,1);
plot(Edc/1000,Vtt/1e9);
xlim([Edcmin Edcmax]/1000);%xlim([-4500 4490]);
ylim([11.4e9 13.5e9]/1e9);
set(gca, 'Position', [0.12 0.75 0.8 0.15],'YTick',[12 13],'YMinorTick','on','ticklength',[0.02 1],'XTick',[-6 -3 0 3 6],'XMinorTick','on')

subplot(2,1,2);
plot(Edc/1000,dgEdc/1e9,'Color',[.5 .5 .5],'LineWidth',1.5);
hold on
plot(Edc(1:1:1200)/1000,ggB0(1:1:1200)/1e9,'Color',[.5 .5 .5],'LineWidth',5);%,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerEdgeColor',[0 0 0],'LineWidth',1);
plot(Edc/1000,dgEdc/1e9+dStark/1e9,'Color',[0 0 0],'LineWidth',1.5);
hold off
legend('A(E),g(E)','Hamiltonian','A(E),g(E),disp(E)','Location','southwest');
ylim([11.185e9 11.207e9]/1e9);
xlim([Edcmin Edcmax]/1000);
set(gca, 'Position', [0.12 0.1 0.8 0.54],'YTick',[11.19 11.20],'YMinorTick','on','ticklength',[0.02 1],'XTick',[-6 -3 0 3 6],'XMinorTick','on')


%% calculate dephasing vs Vt and Edc
% Figure 2b
clearvars;
%addpath Z:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

A=117e6;
B0=0.2*2; %this determine the range of Vt

L=7.5e-9; %half interdot distance
dg=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

% define Hamiltonian matrices %%
% basis =>  nucleus x orbital x electron
% 01) u,e,u
% 02) u,e,d
% 03) u,g,u
% 04) u,g,d
% 05) d,e,u
% 06) d,e,d
% 07) d,g,u
% 08) d,g,d

H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

Edcnoisedim=5;
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=100/sqrt(sum(Edcnoise.^2)/Edcnoisedim);

%  Vt & Edc sweep
dimVt=400;
detmin=60e6;
detmax=350e6;
Vt=sqrt(((28e9+17.2e6-28e6)*B0)^2 + A^2/4)+(detmin:(detmax-detmin)/(dimVt-1):detmax);

dimEdc=300;
Edcmin=-1200;
Edcmax=400;
Edc=Edcmin:(Edcmax-Edcmin)/(dimEdc-1):Edcmax;
Dephasing=zeros(dimVt,dimEdc);
RelRate=zeros(dimVt,dimEdc);

for ii=1:dimVt
    ii
    H_tunel=tensor(identity,Vt(ii)/2*sigma_z,identity);
    for jj=1:dimEdc
        H_Edc=tensor(identity,e*Edc(jj)*L/h*sigma_x,identity);
        
        H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
        
        [H_eVec,H_eVal]=eig(H_0);
        
        nuE=abs(H_eVal(1,1)-H_eVal(3,3)); %these are the eigenvectors that project onto the flipflop close to the ionization point
        
        % sweep Edc noise and average
        Edcnoise=Edc(jj)+((-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax);
        
        for kk=1:Edcnoisedim
            H_Edc=tensor(identity,e*Edcnoise(kk)*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
            [H_eVec,H_eVal]=eig(H_0);
            nuEnoise=abs(H_eVal(1,1)-H_eVal(3,3));
            Dephasing(ii,jj)=Dephasing(ii,jj)+abs(nuE-nuEnoise)/Edcnoisedim;
        end
        
        % Relaxation rate from Palyi's paper
        eo=sqrt(Vt(ii)^2+(e*Edc(jj)*2*L/h)^2);
        pos=-(e*Edc(jj)*2*L/h)/eo;
        eff=sqrt(((28e9*(1+dg*(1+pos)/2)+17.2e6)*B0)^2 + (A*(1-pos)/2)^2);
        gso=A/4*Vt(ii)/eo;
        dso=(eo-eff);
        RelRate(ii,jj)=.49e6/(5.91e9)^3*(eo*(Vt(ii))^2)*(gso/dso)^2;%4*gso^2*Vt(ii)^2*eff^3/eo^3/(eo^2-eff^2)^2;%
    end
end
%%
figure(54675);
subplot(2,2,1)
imagesc(-Edc,Vt*1e-9,log10(Dephasing));colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Vt(GHz)');title('log[Dephasing rate (Hz)]');
set(gca,'Ydir','normal')
subplot(2,2,2)
imagesc(-Edc,Vt*1e-9,log10(RelRate));colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Vt(GHz)');title('log[1/T1 (Hz)]');
set(gca,'Ydir','normal')
subplot(2,2,3)
imagesc(-Edc,Vt*1e-9,log10(Dephasing+RelRate/2));colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Vt(GHz)');title('log[Dephasing + 1/2T1 (Hz)]');
set(gca,'Ydir','normal')

%% Origin
VtO=Vt.'*1e-9;
EdcO=-Edc;
DephasingO=Dephasing;
min(min(DephasingO))
max(max(DephasingO))
RelRateO=RelRate;
min(min(RelRateO))
max(max(RelRateO))

%% Sweep B0 and calculate dephasing and 1/T1 at 2nd order CT
%figure 2e
clearvars;
%addpath Z:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

A=117e6;
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

H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);

B0dim=29;
B0min=0.1;
B0max=0.8;
B0=B0min:(B0max-B0min)/(B0dim-1):B0max;

Edcnoisedim=7;
Edcnoise=-1:2/(Edcnoisedim-1):1;
Edcnoisemax=100/sqrt(sum(Edcnoise.^2)/Edcnoisedim);

dimVt=200;
detmin=160e6;
detmax=400e6;

dimEdc=50;
Edcmin=-570;
Edcmax=-300;
Edc=Edcmin:(Edcmax-Edcmin)/(dimEdc-1):Edcmax;
RelRate=zeros(dimVt,dimEdc);
Deph2CT=zeros(1,B0dim);
Rel2CT=zeros(1,B0dim);

for bb=1:B0dim
    tic;
    bb
    Vt=sqrt(((28e9+17.2e6-28e6)*B0(bb))^2 + A^2/4)+(detmin:(detmax-detmin)/(dimVt-1):detmax);
    H_Znuc=-tensor(17.2e6*B0(bb)/2*sigma_z,identity,identity);
    H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0(bb)/2*sigma_z);
    Dephasing=zeros(dimVt,dimEdc);
    for ii=1:dimVt
        H_tunel=tensor(identity,Vt(ii)/2*sigma_z,identity);
        for jj=1:dimEdc
            H_Edc=tensor(identity,e*Edc(jj)*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
            [H_eVec,H_eVal]=eig(H_0);
            nuE=abs(H_eVal(1,1)-H_eVal(3,3));
            
            % sweep Edc noise and average
            Edcnoise=Edc(jj)+((-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax);
            
            for kk=1:Edcnoisedim
                H_Edc=tensor(identity,e*Edcnoise(kk)*L/h*sigma_x,identity);
                H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
                [H_eVec,H_eVal]=eig(H_0);
                nuEnoise=abs(H_eVal(1,1)-H_eVal(3,3));
                Dephasing(ii,jj)=Dephasing(ii,jj)+abs(nuE-nuEnoise)/Edcnoisedim;
            end
            
            % Relaxation rate from Palyi's paper
            eo=sqrt(Vt(ii)^2+(e*Edc(jj)*2*L/h)^2);
            pos=-(e*Edc(jj)*2*L/h)/eo;
            eff=sqrt(((28e9*(1+dg*(1+pos)/2)+17.2e6)*B0(bb))^2 + (A*(1-pos)/2)^2);
            gso=A/4*Vt(ii)/eo;
            dso=(eo-eff);
            RelRate(ii,jj)=(eo/5.91e9*(Vt(ii)/5.91e9)^2*.49e6)*(gso/dso)^2;%4*gso^2*Vt(ii)^2*eff^3/eo^3/(eo^2-eff^2)^2;%
        end
    end
    [mm,vv]=min(Dephasing);
    [mm,ee]=min(mm);
    Deph2CT(bb)=Dephasing(vv(ee),ee);
    Rel2CT(bb)=RelRate(vv(ee),ee);
    save 1Q-deph_T1_B0_01.mat
    toc;
end
%%
figure(546751);
plot(B0,log10(Deph2CT),B0,log10(Rel2CT/2));
xlabel('B0 (T)');ylabel('log10[Dephasing rate (Hz)]');

%% Origin
BB0=B0';
DDeph2CT=Deph2CT';
DRel2CT2=Rel2CT'/2;

%% flipflop dispersion for three Vt values
% figure 2c
clear all;
addpath Z:\spin-QED\Theory_matlab

clear Edc Vtt ggB0
e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

Edcdim=400;
Edcmin=-400;%332;%
Edcmax=1200;%342;%
Edc=Edcmin:(Edcmax-Edcmin)/(Edcdim-1):Edcmax;
Vt=[11.34e9 11.44e9 11.54e9];

ggB0=zeros(Edcdim,size(Vt,2));

L=7.5e-9; %half interdot distance
A=117e6; %contact hyperfine
B0=0.2*2;
dg=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

orb_e=tensor([1 0],[1 0],[0 1]);
flip_e=tensor([0 1],[0 1],[1 0]);
orbflip_g=tensor([1 0],[0 1],[0 1]);

tic;
for jj=1:size(Vt,2)
H_tunel=tensor(identity,Vt(jj)/2*sigma_z,identity);
for ii=1:Edcdim
    %ii
    H_Edc=-tensor(identity,e*Edc(ii)*L/h*sigma_x,identity);
    H=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
    [H_eVec,H_eVal]=eig(H);
    
    fid1=0;
    fid2=0;
    fid3=0;
    aa=0;
    bb=0;
    cc=0;
    for i=1:8
        AA=abs(orb_e*H_eVec(:,i))^2;
        BB=abs(flip_e*H_eVec(:,i))^2;
        CC=abs(orbflip_g*H_eVec(:,i))^2;
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
    
    ggB0(ii,jj)=abs(H_eVal(bb,bb)-H_eVal(cc,cc));
end
end
toc;

figure(343)
plot(Edc,ggB0/1e9,Edc,(ggB0(:,1)+ggB0(:,2))/2e9,Edc,(ggB0(:,3)+ggB0(:,2))/2e9,Edc,(ggB0(:,3)+ggB0(:,1))/2e9);%,Edc,dgEdc+dStark,Edc,dStark+5.60e9,Edc,dStarkVt+5.60e9
xlim([-400 1200]);%xlim([Edcmin Edcmax]);%
ylim([11.1898 11.195]);
set(gca,'YTick',[11.191 11.192 11.193 11.194],'ticklength',[0.025 1],'XTick',[0 500 1000])%,'XMinorTick','on'
ax = gca;
ax.XAxis.MinorTick = 'on'; ax.XAxis.MinorTickValues = -250:250:1000;
ax.YAxis.MinorTick = 'on'; ax.YAxis.MinorTickValues = 11.1905:0.0005:11.1945;

%% calculate dephasing charge qubit
clear all;
addpath Q:\spin-QED\Theory_matlab\Matlab_functions

e=1.602176565e-19; %C 
h=6.62606957e-34; %J*s

A=117e6;
B0=0.2;

L=7.5e-9; %half interdot distance
dg=-0.002;%relative change in g-factor (delta_g/g)

% Pauli matrices
sigma_x=[0 1;1 0];
sigma_y=[0 -1i;1i 0]; 
sigma_z=[1 0; 0 -1];
identity=eye(2);

H_Hyper=(tensor(sigma_x,identity,sigma_x)+tensor(sigma_y,identity,sigma_y)+tensor(sigma_z,identity,sigma_z))*A/4*tensor(identity,identity/2-sigma_x/2,identity);
H_Znuc=-tensor(17.2e6*B0/2*sigma_z,identity,identity);
H_Zel=tensor(identity,identity+(identity/2+sigma_x/2)*dg,28e9*B0/2*sigma_z);

Edcnoisedim=41;
Edcnoisemax=170;

%  Vt & Edc sweep
dimVt=28;%28*4;
%ii=1:dimVt;
Vtmin=5.9e9;
Vtmax=10e9;
Vt=Vtmin:(Vtmax-Vtmin)/(dimVt-1):Vtmax;
dimEdc=28;%(56+28)*4;
%jj=1:dimEdc;
Edcmin=-50;
Edcmax=30;
Edc=Edcmin:(Edcmax-Edcmin)/(dimEdc-1):Edcmax;
Dephasing=zeros(dimVt,dimEdc);

for ii=1:dimVt
    ii
    H_tunel=tensor(identity,Vt(ii)/2*sigma_z,identity);
    for jj=1:dimEdc
        H_Edc=tensor(identity,e*Edc(jj)*L/h*sigma_x,identity);
        
        H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
        
        [H_eVec,H_eVal]=eig(H_0);
        
        nuE=abs(H_eVal(1,1)-H_eVal(5,5));%nuE=(28e9+17.2e6)*B0*0.9985;%+5e6;
        
        % sweep Edc noise and average
        Edcnoise=Edc(jj)+((-Edcnoisemax):(2*Edcnoisemax)/(Edcnoisedim-1):Edcnoisemax);
        
        for kk=1:Edcnoisedim
            H_Edc=tensor(identity,e*Edcnoise(kk)*L/h*sigma_x,identity);
            H_0=H_Znuc+H_tunel+H_Zel+H_Hyper+H_Edc;
            [H_eVec,H_eVal]=eig(H_0);
            nuEnoise=abs(H_eVal(1,1)-H_eVal(5,5));
            
            Dephasing(ii,jj)=Dephasing(ii,jj)+abs(nuE-nuEnoise)/Edcnoisedim;
        end
    end
end

figure(gcf);
%figure;
imagesc(-Edc,Vt*1e-9,log10(Dephasing));colorbar;%caxis([20e-9 150e-9])
xlabel('Edc(V/m)');ylabel('Vt(GHz)');title('log[Dephasing rate (Hz)]');

