%% simple model %%
% clear all;
addpath \scipnns001.ad.unsw.edu.au\spin-QED\Matlab_functions

muB=9.27400968e-24; %J/T
h=6.62606957e-34; %J*s
hbar=1.054571726e-34; %J*s
e=1.602176565e-19; %C 
e0=8.854187817e-12; %dielectric constant

L=5.5e-9; %half interdot distance

B0=0.5;

% Energy unit is Hz (E/h)
Vt=9.41e9; %double-dot tunel rate
A=108.5e6; %contact hyperfine

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

%only the orbitals are necessary here, can remove the elec and nuc tensors

H_Hyper=(tensor(identity,sigma_x)+tensor(identity,sigma_y)+tensor(identity,sigma_z))*A/4*tensor(identity/2-sigma_x/2,identity);
H_tunel=tensor(Vt/2*sigma_z,identity);
H_Znuc=-tensor(identity,identity);
H_Zeeman=tensor(identity,28e9*B0/2*sigma_z);

Edcmin=-25000;%195e-3;
Edcmax=25000;%220e-3;
Edcdim=1000;
eVal=zeros(Edcdim+1,8);
Hyper=zeros(Edcdim+1,8);
Trans=zeros(Edcdim+1,4);
EEdc=Edcmin:((Edcmax-Edcmin)/Edcdim):Edcmax;
%B0step=0.6e-3;
%BB0=B0min:B0step:B0max;
%B0dim=(B0max-B0min)/B0step+1;
for kk=1:(int16(Edcdim)+1);
    H_Edc=tensor(e*EEdc(kk)*L/h*sigma_x,identity);
    H=H_tunel-H_Edc;
    [H_eVec,H_eVal]=eig(H);
    for ll=1:4;
    eVal(kk,ll)=H_eVal(ll,ll);
    Hyper(kk,ll)=H_eVec(:,1).'*(A*tensor(identity/2-sigma_x/2,identity))*H_eVec(:,1);
    end
    Trans(kk,1)=H_eVal(2,2)-H_eVal(1,1);
%     Trans(kk,2)=H_eVal(6,6)-H_eVal(2,2);
%     Trans(kk,3)=H_eVal(8,8)-H_eVal(3,3);
%     Trans(kk,4)=H_eVal(7,7)-H_eVal(4,4);

    Vtp(kk) = sqrt(Vt^2+(e*2*L*EEdc(kk)/h)^2);
    theta=((e*2*L*EEdc(kk)/h)-sqrt((e*2*L*EEdc(kk)/h)^2+Vt^2))/Vt;
    test(kk)=(theta/sqrt(1+theta^2))^2*A;
end

%% NEMO Hyperfine %%

% bulk wave hyperfine MHz
bulk_hf = 117.6;
bulk_wf_at_donor_site = 0.01774172;
conversion_factor = bulk_hf/bulk_wf_at_donor_site;

% electric_field wf_at_Donor_site in MV/m

data = [2.850    1.7270191684e-02  
        2.900    1.7249676992e-02
        2.950    1.7228366234e-02
        3.000    1.7206274827e-02
        3.050    1.7183334409e-02
        3.100    1.7159504905e-02
        3.150    1.7134716067e-02 
        3.200    1.7108919720e-02
        3.250    1.7082043372e-02
        3.300    1.7054004620e-02
        3.350    1.7024698715e-02
        3.400    1.6994074293e-02
        3.450    1.6961987338e-02
        3.500    1.6928317043e-02
        3.550    1.6892912470e-02
        3.600    1.6855635709e-02
        3.650    1.6816274060e-02
        3.700    1.6774600155e-02
        3.750    1.6730350980e-02
        3.800    1.6683195374e-02
        3.850    1.6632675706e-02
        3.900    1.6578137628e-02
        3.950    1.6518280050e-02
        4.000    1.6449145484e-02
        4.045    1.6355000458e-02
        4.050    1.6336300921e-02
        4.055    1.6312062937e-02
        4.060    1.6278064191e-02
        4.065    1.6224769161e-02
        4.070    1.6126482216e-02
        4.075    1.5893420853e-02
        4.080    1.5038845601e-02
        4.081    1.4606412182e-02
        4.082    1.3956949921e-02
        4.083    1.2967654681e-02
        4.084    1.1487177581e-02
        4.085    9.4726702220e-03
        4.086    7.1710688001e-03
        4.087    5.0752135392e-03
        4.088    3.5006251525e-03
        4.089    2.4368827086e-03
        4.090    1.7434236054e-03
        4.095    4.9808136410e-04
        4.100    2.2160338407e-04
        4.105    1.2297677647e-04
        4.110    7.7493110967e-05
        4.115    5.2996274366e-05
        4.150    1.0255703371e-05
        4.200    2.8995217644e-06
        4.250    1.2539826430e-06
        4.300    6.5956331599e-07
        4.350    3.8849832365e-07
        4.400    2.4639934131e-07
        4.450    1.6464759115e-07
        4.500    1.1436854797e-07
        4.550    8.1861756577e-08
        4.600    6.0006569826e-08
        4.650    4.4850030909e-08
        4.700    3.4067667106e-08
        4.750    2.6233261569e-08
        4.800    2.0438114344e-08
        4.850    1.6085213106e-08
        4.900    1.2771112860e-08
        4.950    1.0219756370e-08
        5.000    8.2355154580e-09
        5.050    6.6778924853e-09
        5.100    5.4454025401e-09
        5.150    4.4631846684e-09
        5.200    3.6752784165e-09
        5.250    3.0393804567e-09
        5.300    2.5234087592e-09
        5.350    2.1026546962e-09
        5.400    1.7579611789e-09
        5.450    1.4743364973e-09       
    ];

field_MV_m = data(:,1);
wf_at_donor_site = data(:,2);
hyperfine_MHZ = wf_at_donor_site*conversion_factor; 

hyperfine_field_gradient = diff(hyperfine_MHZ)./diff(field_MV_m);

%% NEMO energies %%

% bulk wave hyperfine MHz
bulk_hf = 117.6;
bulk_wf_at_donor_site = 0.01774172;
conversion_factor = bulk_hf/bulk_wf_at_donor_site;

% electric_field E0 E1 in MV/m



data = [2.850    1.128759968514e+00  1.140271520596e+00      
        2.900    1.129507206752e+00  1.141006377628e+00    
        2.950    1.130254147620e+00  1.141734207897e+00    
        3.000    1.131000727798e+00  1.142872319290e+00   
        3.050    1.131746947380e+00  1.142990982953e+00 
        3.100    1.132492776036e+00  1.143265612659e+00       
        3.150    1.133238232718e+00  1.143484716093e+00         
        3.200    1.133983266652e+00  1.143693700286e+00       
        3.250    1.134727890722e+00  1.143898354175e+00       
        3.300    1.135472077363e+00  1.144100266525e+00       
        3.350    1.136215828046e+00  1.144300072680e+00     
        3.400    1.136959078977e+00  1.144498100027e+00        
        3.450    1.137701830529e+00  1.144694539238e+00    
        3.500    1.138444080679e+00  1.144889512110e+00    
        3.550    1.139185770336e+00  1.145083110285e+00    
        3.600    1.139926863999e+00  1.145275403711e+00   
        3.650    1.140667358188e+00  1.145466450225e+00    
        3.700    1.141407193478e+00  1.145656296363e+00    
        3.750    1.142146318485e+00  1.145844985930e+00        
        3.800    1.142884679709e+00  1.146032557591e+00    
        3.850    1.143622252689e+00  1.146219041829e+00    
        3.900    1.144358934956e+00  1.146404482536e+00    
        3.950    1.145094624854e+00  1.146588919446e+00    
        4.000    1.145829150808e+00  1.146772431529e+00   
        4.045    1.146488855836e+00  1.146937055890e+00   
        4.050    1.146561995291e+00  1.146955366146e+00
        4.055    1.146635093700e+00  1.146973702910e+00    
        4.060    1.146708115645e+00  1.146992093695e+00    
        4.065    1.146780988034e+00  1.147010570126e+00   
        4.070    1.146853644954e+00  1.147029234349e+00   
        4.075    1.146925812015e+00  1.147048339182e+00    
        4.080    1.146996463022e+00  1.147068892307e+00 
        4.081    1.147010115755e+00  1.147073485921e+00
        4.082    1.147023395922e+00  1.147078432002e+00
        4.083    1.147036125683e+00  1.147083912953e+00
        4.084    1.147048057384e+00  1.147090218205e+00
        4.085    1.147058792435e+00  1.147097692420e+00
        4.086    1.147068053390e+00  1.147106642564e+00
        4.087    1.147075801994e+00  1.147117094960e+00
        4.088    1.147082305316e+00  1.147128789632e+00
        4.089    1.147087924589e+00  1.147141372855e+00
        4.090    1.147092956886e+00  1.147154527470e+00    
        4.095    1.147114306603e+00  1.147223993731e+00    
        4.100    1.147133560900e+00  1.147295121740e+00
        4.105    1.147152231571e+00  1.147365807788e+00   
        4.110    1.147170662098e+00  1.147432555575e+00    
        4.115    1.147188963464e+00  1.147479628830e+00    
        4.150    1.147315831236e+00  1.147620305001e+00 
        4.200    1.147495639821e+00  1.147805006541e+00   
        4.250    1.147674419525e+00  1.147988340759e+00    
        4.300    1.147852282409e+00  1.148170701730e+00    
        4.350    1.148029272176e+00  1.148352168780e+00    
        4.400    1.148205413850e+00  1.148532778740e+00    
        4.450    1.148380731193e+00  1.148712554604e+00      
        4.500    1.148555237509e+00  1.148891513942e+00   
        4.550    1.148728953277e+00  1.149069678972e+00  
        4.600    1.148901894166e+00  1.149247066141e+00   
        4.650    1.149074077197e+00  1.149423689385e+00    
        4.700    1.149245505859e+00  1.149599556043e+00   
        4.750    1.149416204857e+00  1.149774689104e+00   
        4.800    1.149586186157e+00  1.149949103538e+00           
    ];

field_MV_mm = data(:,1);
E_0 = data(:,2);
E_1 = data(:,3);

CB = 1.085716951844e+00 + 45.6e-3;
norm_conduction_band = CB + ((15.20666e-9)*field_MV_mm*(1e6));
E_0 = E_0 - norm_conduction_band;
E_1 = E_1 - norm_conduction_band;

%% plot all %%

figure(gcf);
subplot (2,1,1)
H1=plot(field_MV_m, hyperfine_MHZ,EEdc/1e6+4.0856,Hyper/1e6);
hold on;
plot(EEdc/1e6+4.0856, test/1e6, 'Color', 'r');
%hold off
xlabel('Electric Field Ez (MV/m)');
ylabel('Hyperfine Splitting (MHz)');
axis([4.062 4.111 -3 115])
set(H1(1),'Marker','o','LineStyle','none','LineWidth',2,'Color',[1 0 0]);
set(H1(2),'Marker','none','LineStyle','-','LineWidth',2,'Color',[0 1 1]);
% subplot(4,1,2);
% plot(EEdc/1e6+4.085,Hyper/1e6);
%axis([Edcmin Edcmax 0 110e6])
%set(gca,'Xlim',[Edcmin Edcmax]);
subplot(2,1,2)
H2=plot(field_MV_mm, (E_1 - E_0)*2.418e14,EEdc/1e6+4.0856,Vtp);
xlabel('Electric Field Ez  (MV/m)');
ylabel('Energy separation between ground and excited states (Hz)');
xlim([4.047 4.125])
set(H2(1),'Marker','o','LineStyle','none','LineWidth',2,'Color',[1 0 0]);
set(H2(2),'Marker','none','LineStyle','-','LineWidth',2,'Color',[0 0 1]);
% subplot(3,1,3);
% plot(EEdc,Trans);figure(gcf);
% set(gca,'Xlim',[Edcmin Edcmax]);


% subplot(4,1,1);
% plot(EEdc,eVal);figure(gcf);
% set(gca,'Xlim',[Edcmin Edcmax]);
% set(gca,'Ylim',[7.9e8 9.1e8],'Xlim',[0.1989 0.2051])
% subplot(4,1,2);
% imagesc(abs(H_eVec).^2);colorbar;figure(gcf);
% set(gca,'Xlim',[2.5 8.5])

%%
h=figure(11)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
c=colormap(parula)
H2=plot(field_MV_mm, (E_1 - E_0)*2.418e14/1e9,EEdc/1e6+4.0856,Vtp/1e9);
xlabel('Electric field E_z  (MV/m)', 'fontsize', 18);
ylabel('\epsilon_o (GHz)');
xlim([4.061 4.11])
set(H2(1),'Marker','o','LineStyle','none','LineWidth',2,'Color',c(1, :));
set(H2(2),'Marker','none','LineStyle','-','LineWidth',2,'Color',c(25,:));
set(gca,'FontSize',16, 'LineWidth', 3);
%legend('Nemo3D similation', 'theory')
set(gca,'color','none')
print('chargequbit_dispersion.pdf','-r0','-dpdf')
%saveas(gcf,'chargequbit_dispersion.pdf')