% EQUIVALENT SCRIPT FOR PYTHON COMPARISON
% same settings as 2LAY_BASE_CODE.PY

% std normal size distrib with V_eff = 0.2 vs V_eff = 0.6 
% static ci 2.5e6, nshell = 1.10, ncore=1.02

clear

savedir='/Users/user/Documents/WORK/HYDRO/'
int_val=56; % 56 for 400:5:900; 276 for 400:1:900
filename = 'test_for_py_comparison_effvar2' % name for resulting IOP mat:
lambda=[0.400:0.005:0.900]';
D_eff= [2,8,16]';
V_eff=2; % Effective variance
Vs=0.2;	% Shell as chloroplast volume   
ci=2.5e6

%*********************************************************************************
%		Load and set up refractive index  and cell geometry related data
%*********************************************************************************
cd /Users/user/Documents/WORK/2LAYER/2lay/
%Get name of input file for lambda,n, k and PSD
load 501nm_extended_e1701000

Vc=1-Vs;     %Core as cytoplasm volume
%Calculate relative radius based on outer chloroplast volume
FR=(1-Vs).^(1/3);

%Define real refractive index of surrounding medium and calculate effects on wavelength
nmedia=1.334;
wavelength=lambda./nmedia; %NB wavelength now in MICRONS - no it's not??
wvno=2.*pi./wavelength;
   
%Load and calculate refractive index of spheres(input file has relative n)
kcore=interp1(RIs(:,6),RIs(:,1),lambda,'linear'); kcore(find(kcore<1e-15))=1e-15;
kshell_base=interp1(RIs(:,6),RIs(:,3),lambda,'linear');

%Normalise n'(675) to 0.019942, based on ci=2.5e6 kg m-3 and chloroplast volume=0.2   
kshell_norm=(6.75e-7./nmedia).*(0.027.*ci./Vs)./(4.*pi);
   
%find kshell val at lambda 675: 58 for 390 to 900 nm, normalise to that.
%286 for 511 nm. 276 for 501 nm.(400:1:900), 56 for 400:900 nm at 5nm.
kshell=kshell_base.*(kshell_norm./kshell_base(int_val)); 
kshell(find(kshell<1e-15))=1e-15;
   
nshell=1.10+imag(hilbert(kshell));% real part of refractive index, default is 1.10, can be varied from 1.08 - 1.14
ncore=1.02+imag(hilbert(kcore)); % default is 1.02

khom=kcore.*Vc+kshell.*Vs;
nhom=ncore.*Vc+nshell.*Vs;
   
%Define inital estimate of k and n for shell
mshell=nshell-kshell.*i;
mcore=ncore-kcore.*i;
mhom=nhom-khom.*i;

%*********************************************************************************
%	Set up inital size distribution data
%*********************************************************************************

%Set up size vector for radius and alpha   
psd(:,1)=[1:1:100]';
deltad=1; deltadm=1;
      
%*********************************************************************************
%	Set up angular specifications for Dmilay scattering routine
%*********************************************************************************
   
%Set up angles for VSF determination, backscattering limits, and delta(theta).
theta901=[0:0.1:90]';
nang=901;
angles=cos(deg2rad(theta901));
%Now set up full angle vector to pi for integration purposes
theta=cat(1,deg2rad(theta901),flipud(pi-deg2rad(theta901(1:900))));
dtheta=diff(theta); dtheta(1801)=theta(1801)-theta(1800);
back1=find(theta==(pi/2));
back2=find(theta==(pi)); 
   
%*********************************************************************************
%	Run Dmilay for all sizes and wavelengths
%*********************************************************************************
cd /Users/user/Documents/WORK/2LAYER/2lay/  
%Start loop to cycle through wavelengths
[ni nj]=size(kcore); 

	for nii=1:ni;
        nii

	    %start loop for summation through alpha vector for beta caculation
		for jj=1:length(psd);
      
        %***************************************************************************
		%	Calculate Aden-Kerker two layered efficiency factors 
		%***************************************************************************		

		%Set up for calculation of beta using DiMilay mex call
		[Qcro(jj),Qbro(jj),g,Qbs,m1,m2,s21,d21]=Lisl_mex_clean((psd(jj,1).*FR)./2,psd(jj,1)./2,wvno(nii),mshell(nii),mcore(nii),angles,901,901);
		%Assume m1 and m2 are i1 and i2, calculate angular intensity parameter i (Morel & Bricaud Can.Bull. Fish.Aq.Sci, 1986)
		%First need to create appropriate intensity vectors due to Dimilay strange angular returns
		M1=[m1(:,1);flipud(m1(1:900,2))]; M2=[m2(:,1);flipud(m2(1:900,2))];
		II(jj,:)=(M1+M2)'./2;

		%Calculate raw phase function [M&B '86, eq. 18] NB this does not produce
		%phase function that satisfies normalisation conditions
        alpha(jj)=2.*pi*(psd(jj,1)./2)./wavelength(nii);
		phaseMB(jj,:)=II(jj,:)./(pi.*Qbro(jj).*alpha(jj).^2);
		checkMB(jj,:)=2*pi.*sum((phaseMB(jj,:)'.*sin(theta)).*dtheta);

		%Calculate backscattering probability from normalised phase function
		bbprob(jj)=(2*pi.*sum((phaseMB(jj,(back1:back2))'.*sin(theta(back1:back2)).*dtheta(back1:back2))));
		%Calculate backscattering efficiency factor for full size range using beta
		Qbbro(jj)=Qbro(jj).*bbprob(jj);
        clear S1 S2 m1 m2 s21 d21 M1 
     
       
        end; %End of PSD loop
      
       
d_alpha(nii)=diff(alpha(1:2));

 %*********************************************************************************
 %	Loop through Deff determined Standard size distributions of PSD integration
 %*********************************************************************************
	
        %Set up loop to cycle through sizes
        for jjj=1:length(D_eff);
         
        %Calculate Standard distribution based on effective diameter   
        psd(:,2)=1e20.*(psd(:,1)./2).^((1-3.*V_eff)./V_eff).*exp((-psd(:,1)./2)./((D_eff(jjj)./2).*V_eff));
   
        %Generate meter based size distribution for integrating purposes
        psdm(:,1)=psd(:,1)./1e6;
        psdm(:,2)=psd(:,2).*1e3; 
        civol=pi./6.*sum(psdm(:,2).*psdm(:,1).^3.*deltadm);
        psdm(:,2)=psdm(:,2).*(1./(civol.*ci));
        psdvol=pi./6.*sum(psdm(:,2).*psdm(:,1).^3.*deltadm);

        %Find index of D_eff in psd vector
        D_eff_i=find(psd(:,1)==D_eff(jjj));
   
        %*****************************************************************************
        %	Size distribution integrations for efficiency factors and IOPS
        %*****************************************************************************
   
        %determination of PSD integral for two layered geometry
        Qc(jjj,nii)=sum(Qcro'.*psdm(:,2).*(psdm(:,1).^2).*deltadm)./sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);
        Sigma_c(jjj,nii)=pi./4.*Qc(jjj,nii).*sum((psdm(:,1).^2).*deltadm);
        c(jjj,nii)=pi/4.*Qc(jjj,nii).*sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);
        Qb(jjj,nii)=sum(Qbro'.*psdm(:,2).*(psdm(:,1).^2).*deltadm)./sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);
        Sigma_b(jjj,nii)=pi./4.*Qb(jjj,nii).*sum((psdm(:,1).^2).*deltadm);
        b(jjj,nii)=pi/4.*Qb(jjj,nii).*sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);
        Qbb(jjj,nii)=sum(Qbbro'.*psdm(:,2).*(psdm(:,1).^2).*deltadm)./sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);
        Sigma_bb(jjj,nii)=pi./4.*Qbb(jjj,nii).*sum((psdm(:,1).^2).*deltadm);
        bb(jjj,nii)=pi/4.*Qbb(jjj,nii).*sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);
        Qa(jjj,nii)=Qc(jjj,nii)-Qb(jjj,nii);
        a(jjj,nii)=c(jjj,nii)-b(jjj,nii);
        Sigma_a(jjj,nii)=pi./4.*Qa(jjj,nii).*sum((psdm(:,1).^2).*deltadm);
        bbtilde(jjj,nii)=bb(jjj,nii)./b(jjj,nii);
   
        %Determination of PSD integral for beta, VSF and integration check
        %(Morel & Bricaud Can.Bull. Fish.Aq.Sci, 1986).
            for ai=1:(nang*2-1);
		% betabar(ai)=sum(phaseMB(:,ai).*psdm(:,2).*(psdm(:,1).^2).*deltadm)./sum(psdm(:,2).*(psdm(:,1).^2).*deltadm);     
        % not using the above any more, was a problem (typo?) - from euqn perhaps based on euqn 17 in M&B
            betabar(ai)=(1./pi).*(sum(II(:,ai).*psdm(:,2).*d_alpha(nii))./sum(Qbro'.*psdm(:,2).*alpha.^(2)'.*d_alpha(nii)));
            VSF_1(ai)=betabar(ai).*b(jjj,nii);   
            end
            
        checkbar=(2*pi.*sum(betabar'.*sin(theta).*dtheta));
        PF(jjj,nii,:)=betabar;
        VSF(jjj,nii,:)=VSF_1;
        PF_check(jjj,nii,:)=checkbar;
        b_check(jjj,nii)=2*pi.*sum((squeeze(VSF(jjj,nii,:)).*sin(theta)).*dtheta);
        bb_check(jjj,nii)=(2*pi.*sum((squeeze(VSF(jjj,nii,(back1:back2))).*sin(theta(back1:back2)).*dtheta(back1:back2))));
          
        %**************************************************************************************************
        %	Package effect calculations
        %**************************************************************************************************
      
        %Set up all phase lag and optical thickness parameters
        %Calculate acm (absorption celular material) for all layers
        acm_core=(4.*pi.*kcore(nii))./(wavelength(nii).*1e-6);
        acm_shell=(4.*pi.*kshell(nii))./(wavelength(nii).*1e-6);
        acm_hom=(4.*pi.*khom(nii))./(wavelength(nii).*1e-6);
        q=(D_eff(jjj)./2.*FR)./(D_eff(jjj)./2); %or is this Deff divided by 2??

        %Qas2 - package effect formulation with 2 layered Aden_Kerker derived Qa   
        %Qas2(jjj,nii)=3./2.*(Qa(jjj,nii)./((acm_core.*q.^3+acm_shell.*(1-q.^3)).*2e-6.*D_eff(jjj)));
        Qas21(jjj,nii)=3./2.*(Qa(jjj,nii)./((acm_core.*q.^3+acm_shell.*(1-q.^3)).*1e-6.*D_eff(jjj)));
        %Qas1 - package effect formulation with homogenous Mie derived Qa   
      
        %Direct volume equivalenwikipediat determination of package effect
        a_sol(jjj,nii)=psdvol.*Vc.*acm_core+psdvol.*Vs.*acm_shell;
        a_solm(jjj,nii)=psdvol.*acm_hom;
        Qastar2_dir(jjj,nii)=a(jjj,nii)./a_sol(jjj,nii);
        % Hayley uses this one Qastar2_dir (for input to fluorescence term in
        % algorithm)
      
        end; %End of effective radius loop
   
    end; %End of wavelength loop

eval(['cd ' savedir])
save(filename)    
