module DielectricFunctions

using Dierckx
using DelimitedFiles
using LinearAlgebra
using InitializeSystem


export drude
export Scorr
export Diel,DrudeAu
export alphap
export alphap1
export uspec
export PlotS
export Flm
export lte
export LtzSiC

export alphaE
export alphaB

export ReadExp, GetDrudeParameters
export h, hbar, c, wpAu, gammaAu,hbar_eV,h_eV,kB_eV
export wpAg, gammaAg, vf
export path1,kB, PlotDielectric
export wpAl,wpCu,gammaAl,gammaCu
export dtr,rtd,KD,lte,drude2

#Plank constant [h] = ev/s
  const h_eV    = 4.135667e-15
  #plank constant hbar
  const hbar_eV = h_eV/(2*pi)
  #speed of light in nm/sec
  const c       = 2.997e17
  #Au Drude plasma frequency in eV
  const  wpAu   = 8.55
  #AL Drude plasma frequency in eV
  const  wpAl   = 15.8
  #Cu Drude plasma frequency in eV
  const  wpCu   = 10.8
  #Au damping constant in eV (gamma*hbar)
  const gammaAu = wpAu*0.0126
  #Ag Drude plasma frequency in eV
  const wpAg    = 9.20
  #Ag damping constant in eV (gamma*hbar)
  const gammaAg = wpAg*0.00188
  #Al damping constant in eV (gamma*hbar)
  const gammaAl = wpAl*0.04
  #Cu damping constant in eV (gamma*hbar)
  const gammaCu = wpCu*0.00225
  #Fermi velocity in nm/sec
  const vf      = 1.4e15
  #InputParameters
  path1=pwd()
  #Boltzman [KBoltz] = ev/K
  const kB_eV   = 8.6173303e-5
  const h       = 6.62607015e-34
  const hbar    = 1.054571800e-34
  const kB      = 1.3806485e-23

# --------------------------------
#       InSb
# --------------------------------

qe = 1.6e-19;          # C, carga del electron
me = 9.10938188e-31;   # kg, masa del electron
e0 = 8.854e-12;        # F/m

# ------ InSb -----------
N_insb     =  1e26;                                       # m^-3
m_eff_insb = 0.0153 * me;
mu_insb    = 18.03 / (1 + (N_insb/1e6/3e17)^0.68);      # m^2/V/s

omega_insb = sqrt(N_insb*qe*qe/e0/m_eff_insb);            # rad/s
t_insb     = m_eff_insb * mu_insb / qe;
e_inf      = 15.68;
lto(lambda) = lte(lambda)/hbar_eV 
eps_insb(lambda)   = e_inf - omega_insb^2 / (lto(lambda) + 1*im/t_insb) / lto(lambda);

global wpInSb    = hbar_eV*omega_insb
global gammaInSb = hbar_eV/t_insb
global EinfInSb  = e_inf

##########################################
#### Converts from wavelength(nm) to eV ##
##########################################
KD(i,j)     = ifelse(i == j ,1,0)
dtr(degree) = degree*(pi/180)
rtd(radian) = radian*(180/pi)
lte(l)      = h_eV*c/l
######################
## Read ExpData      ##
#######################

function ReadExp(nk,i)

  if(nk==1 && i==1)

  intemp=string(path1,"/InputData/jc_Aun2.txt")

  filen=open(intemp)

  return readdlm(filen)

  close(filen)

  elseif(nk==1 && i==2)

  intemp=string(path1,"/InputData/jc_Auk2.txt")
  filek=open(intemp)

  return readdlm(filek)
  close(filek)


  elseif(nk==2 && i==1)

  intemp=string(path1,"/InputData/jc_Agn2.txt")
  filen=open(intemp)

  return readdlm(filen)

  close(filen)

  elseif(nk==2 && i==2)

  intemp=string(path1,"/InputData/jc_Agk2.txt")

  filek=open(intemp)

  return readdlm(filek)
  close(filek)

elseif(nk==3 && i==1)

  intemp=string(path1,"/InputData/Al_n.txt")
  filek=open(intemp)

  return readdlm(filek)
  close(filek)

elseif(nk==3 && i==2)

  intemp=string(path1,"/InputData/Al_k.txt")
  filek=open(intemp)

  return readdlm(filek)
  close(filek)
  elseif(nk==4 && i==1)

  intemp=string(path1,"/InputData/Cu_n.txt")
  filek=open(intemp)

  return readdlm(filek)
  close(filek)

elseif(nk==4 && i==2)

  intemp=string(path1,"/InputData/Cu_k.txt")
  filek=open(intemp)

  return readdlm(filek)
  close(filek)

  end

end


########################################
## Drude Size correction constant
########################################
gsize(a) = vf*hbar/a


function DielParameters(Specie )
  if Specie == "Au" || Specie=="AuExp"
    return (gammaAu,wpAu)
  elseif Specie == "Ag" || Specie=="AgExp" || Specie == "AgPalik"
    return (gammaAg,wpAg)
  elseif Specie == "Cu" || Specie=="CuExp"
    return (gammaCu,wpCu)
  elseif Specie == "Al" || Specie=="AlExp"
    return (gammaAl,wpAl)
  end

end


##################################
##################################
function DrudeAu(omega_eV)
    EInf    = 5.9752
    wpd     = 8.8667
    gammad  = 0.03799
    s1      = 1.76
    wp1l    = 3.6
    gamma1l = 1.3
    s2      = 0.952
    wp2l    = 2.8
    gamma2l = 0.737
    
    term1 = wpd^2/(omega_eV^2 + im*gammad*omega_eV)
    term2 = s1*wp1l^2/((omega_eV^2 - wp1l^2)+im*gamma1l*omega_eV)
    term3 = s2*wp2l^2/((omega_eV^2 - wp2l^2)+im*gamma2l*omega_eV)
    
  return EInf - term1 - term2 - term3

end


###############################
## Lorentzian
###############################
function LtzSiC(lambda)
    omegaL = 969
    omegaT = 793
    gamma  = 4.76
    Einf   = 6.7
    lambdacm = lambda*1.0e-7
 
    wnumber= 1.0/ lambdacm
    
    return Einf*(1.0 + (omegaL^2 - omegaT^2 ) /(omegaT^2-wnumber^2 - im*gamma*wnumber) )
end



function LtzSiO2(lambda)
    
  lto(lambda)     = 2pi*c / lambda
    
  w               = lto(lambda)  
  A1              = 8.2736e+13;
  w01             = 8.54484e+13;
  G1              = 8.46448e+12;
  A2              = 1.58004e+14;
  w02             = 2.029e+14;
  G2              = 1.06449e+13;
  A3              = 3.39786e+13;
  w03             = 1.51198e+14;
  G3              = 8.33205e+12;
  EpsInf          = 2.03843;

  return EpsInf + A1*A1/(w01*w01 - w*w - im*w*G1) + A2*A2/(w02*w02 - w*w - im*w*G2) + A3*A3/(w03*w03 - w*w - im*w*G3)

end

###################################
## Drude and size correction are evaluated in eV
## Drude Function
##################################
function drude2(lambda,wp,gdrude)

  return 1-wp^2/(lte(lambda)*(lte(lambda)+im*gdrude))

end

function drude(lambda,Specie)

  gdrude=DielParameters(Specie)[1]
  wp= DielParameters(Specie)[2]
  return 1-wp^2/(lte(lambda)*(lte(lambda)+im*gdrude))

end

###################################
## Size Correction
##################################
function Scorr(lambda,Specie, a )

  gdrude = DielParameters(Specie)[1]
  wp = DielParameters(Specie)[2]
  return  1-wp^2/(lte(lambda)*(lte(lambda)+im*gdrude+im*gsize(a)))

end

##################################
## Creating Dielectric Function ##
## and Spectral Variable        ##
## Experimental function evaluated in lambdas nm
##################################


nexpdatAu=ReadExp(1,1)
kexpdatAu=ReadExp(1,2)
int1Au=Spline1D(nexpdatAu[:,1],nexpdatAu[:,2])
int2Au=Spline1D(kexpdatAu[:,1],kexpdatAu[:,2])

nexpdatAg=ReadExp(2,1)
kexpdatAg=ReadExp(2,2)
int1Ag=Spline1D(nexpdatAg[:,1],nexpdatAg[:,2])
int2Ag=Spline1D(kexpdatAg[:,1],kexpdatAg[:,2])

expdatAgPalik = readdlm(string(path1,"/InputData/Palik_Ag.txt"),  skipstart = 2)
int1AgPalik=Spline1D(expdatAgPalik[:,1].*1e3,expdatAgPalik[:,2])
int2AgPalik=Spline1D(expdatAgPalik[:,1].*1e3,expdatAgPalik[:,3])

nexpdatAl=ReadExp(3,1)
kexpdatAl=ReadExp(3,2)
int1Al=Spline1D(nexpdatAl[:,1],nexpdatAl[:,2])
int2Al=Spline1D(kexpdatAl[:,1],kexpdatAl[:,2])

nexpdatCu=ReadExp(4,1)
kexpdatCu=ReadExp(4,2)
int1Cu=Spline1D(nexpdatCu[:,1],nexpdatCu[:,2])
int2Al=Spline1D(kexpdatCu[:,1],kexpdatCu[:,2])

# This data is wavelength in microns, real(eps) imag(eps)
dataAg_Ana = readdlm(string(path1,"/InputData/Ag_jc_ana.tab"),  skipstart = 1)

intAg_ana_rdiel  = Spline1D(dataAg_Ana[:,1]*1E3,dataAg_Ana[:,5])
intAg_ana_idiel  = Spline1D(dataAg_Ana[:,1]*1E3,dataAg_Ana[:,6])

######### in this case [x] = nm is wavelenght #########
function Diel(lambda_nm ,specie,a)
  a = a
    
if specie == "Au"
   return (evaluate(int1Au,lambda_nm))^2-(evaluate(int2Au,lambda_nm))^2+im*(2*evaluate(int2Au,lambda_nm)*evaluate(int1Au,lambda_nm)-drude(lambda_nm,specie)+Scorr(lambda_nm,specie,a) )
elseif specie == "AuDrude"
   return  drude(lambda_nm,"Au")   
 elseif specie == "AuExp"
   return (evaluate(int1Au,lambda_nm))^2-(evaluate(int2Au,lambda_nm))^2+im*(2*evaluate(int2Au,lambda_nm)*evaluate(int1Au,lambda_nm))
 elseif specie == "Ag"
   return (evaluate(int1Ag,lambda_nm))^2-(evaluate(int2Ag,lambda_nm))^2+im*(2*evaluate(int2Ag,lambda_nm)*evaluate(int1Ag,lambda_nm))-drude(lambda_nm,specie)+Scorr(lambda_nm,specie,a)
 elseif specie == "AgExp"
   return (evaluate(int1Ag,lambda_nm))^2-(evaluate(int2Ag,lambda_nm))^2+im*(2*evaluate(int2Ag,lambda_nm)*evaluate(int1Ag,lambda_nm))
 elseif specie == "AgPalik"
       return (evaluate(int1AgPalik,lambda_nm))^2-(evaluate(int2AgPalik,lambda_nm))^2+im*(2*evaluate(int2AgPalik,lambda_nm)*evaluate(int1AgPalik,lambda_nm))   
elseif specie == "AgDrude"
   return  drude(lambda_nm,"Ag") 
elseif specie == "AgAna"
        return evaluate(intAg_ana_rdiel,lambda_nm) + im * evaluate(intAg_ana_idiel,lambda_nm)
elseif specie == "Al"
   return (evaluate(int1Al,lambda_nm))^2-(evaluate(int2Al,lambda_nm))^2+im*(2*evaluate(int2Al,lambda_nm)*evaluate(int1Al,lambda_nm))-drude(lambda_nm,specie)+Scorr(lambda_nm,specie,a)
 elseif specie == "AlExp"
   return (evaluate(int1Al,lambda_nm))^2-(evaluate(int2Al,lambda_nm))^2+im*(2*evaluate(int2Al,lambda_nm)*evaluate(int1Al,lambda_nm))
 elseif specie == "Cu"
  return (evaluate(int1Cu,lambda_nm))^2-(evaluate(int2Cu,lambda_nm))^2+im*(2*evaluate(int2Cu,lambda_nm)*evaluate(int1Cu,lambda_nm))-drude(lambda_nm,specie)+Scorr(lambda_nm,specie,a)
elseif specie == "CuExp"
  return (evaluate(int1Cu,lambda_nm))^2-(evaluate(int2Cu,lambda_nm))^2+im*(2*evaluate(int2Cu,lambda_nm)*evaluate(int1Cu,lambda_nm))
elseif specie == "SiC"
  return LtzSiC(lambda_nm)
elseif specie == "SiO2"
  return LtzSiO2(lambda_nm)
elseif specie == "AuDrude2"
  return DrudeAu(lte(lambda_nm))        
elseif specie == "InSb"
  return eps_insb(lambda_nm)
end
    
end


function uspec(lambda,specie,Emed,a)
   return 1/(1-Diel(lambda,specie,a)/Emed)
end

function uspec(lambda,wp,gdrude,Emed,a)
   return 1/(1-drude2(lambda,wp,gdrude)/Emed)
end

N0(l) =  l/(2*l+1);

function alphap(lambda,l,specie,Emed,a)
    return N0(l)*a^(2*l+1)/(N0(l)-uspec(lambda,specie,Emed,a) )

end

function PlotDielectric(Specie,a)
    lambda_list = range(200,800, length = 100)
    return plot(lambda_list, hcat( real.(Diel.(lambda_list ,Specie,a)), imag.(Diel.(lambda_list ,Specie,a))), lw  = 3,
lab = ["Real" "Imag"],title = "$Specie Dielectric Function",xlabel = "wavelength (nm)")
end


#######################
### Taken from
### Radiate Heat transfer in many body systems: 
### coupled electric and coupled magnetic approach
####################
Jsp1(x)    = sin(x)/x^2 - cos(x)/x
Ysp1(x)    = -cos(x)/x^2 - sin(x)/x 
Hsp1(x)    = Jsp1(x) + im*Ysp1(x)
DJsp1(x)   = 2*cos(x)/x^2 - 2*sin(x)/x^3 + sin(x)/x
DYsp1(x)   = 2*cos(x)/x^3 - cos(x)/x + (2*sin(x))/x^2
DHsp1(x)   = DJsp1(x)  + im*DYsp1(x)


function alphaE_QS(omega_eV,specie,Emed,a)
     return  a^3* (Diel(lte(omega_eV),specie,a)- Emed) / (Diel(lte(omega_eV),specie,a)+ 2*Emed)

end

#######################
### Taken from
### Radiate Heat transfer in many body systems: 
### coupled electric and coupled magnetic approach
####################
function alphaE(omega_eV,specie,Emed,a)
    
   # Check Bohren hofmann Eq. 4.53
    k      = 2pi*sqrt(Emed)/lte(omega_eV)
    x      = k*a
    m      = sqrt(Diel(lte(omega_eV),specie,a)/Emed)
    y      = m*x
    coef1  = (3*im/(2*(k)^3))
    
    Der1   = x*DJsp1(x) + Jsp1(x)
    Der2   = y*DJsp1(y) + Jsp1(y)
    Der3   = x*DHsp1(x) + Hsp1(x)
    Der4   = Der2


    a1_num = m^2*Jsp1(y)*Der1 - Jsp1(x)*Der2
    a1_den = m^2*Jsp1(y)*Der3 - Hsp1(x)*Der4

    
    return coef1*(a1_num/a1_den)

end


function alphaB(omega_eV,specie,Emed,a)
    
    k      = 2pi*sqrt(Emed)/lte(omega_eV)
    x      = k*a*sqrt(Emed)
    m      = sqrt(Diel(lte(omega_eV),specie,a))/sqrt(Emed)
    y      = m*x
    coef1  = (3*im/(2*(k)^3))
    
    Der1   = x*DJsp1(x)  + Jsp1(x)
    Der2   = y*DJsp1(y) + Jsp1(y)
    Der3   = x*DHsp1(x) + Hsp1(x)
    Der4   = Der2


    a1_num = Jsp1(y)*Der1 - Jsp1(x)*Der2
    a1_den = Jsp1(y)*Der3 - Hsp1(x)*Der4
    
    return coef1*(a1_num/a1_den)

end

function GetDrudeParameters(SP::SystemParameters)
    if SP.Specie == "Ag" || SP.Specie == "AgExp" || SP.Specie == "AgDrude"
        return (wpAg,gammaAg,1)
    elseif SP.Specie == "Au" || SP.Specie == "AuExp" || SP.Specie == "AuDrude"
        return (wpAu,gammaAu,1)
    elseif SP.Specie == "Al" || SP.Specie == "AlExp" || SP.Specie == "AlDrude"
        return (wpAl,gammaAl,1)
    elseif SP.Specie == "InSb"
        return (wpInSb,gammaInSb,EinfInSb)
    end
end



end #endModule
