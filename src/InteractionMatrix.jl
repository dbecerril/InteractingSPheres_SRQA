module InteractionMatrix

using DielectricFunctions
using InitializeSystem
using GSL
using LinearAlgebra

export CartToSphere
export Aelement, Aelement_periodic
export Helement
export Hmat
export KD,NO,SepVec
export Hmat_periodic,Hmat_PiPj
export Hmat_periodic_temp

function CartToSphere(Rs)
    R = norm(Rs)
    x = Rs[1]
    y = Rs[2]
    z = Rs[3]

#### defining the origin ######
    if R == 0.0
        phi   = 0.0
        theta = 0.0
    end

####### Defining theta
#    if x==0.0 && y==0.0
#        if z > 0.0
#            theta = 0.0
#        elseif z < 0.0
#            theta = pi
#        end
#    elseif z>0
        theta = acos(z/R)
#    else 
#        theta = acos(-z/R)
#    end

###### Defining phi
#    if x==0.0 && y>0.0
#            phi=  0.5*pi
#    elseif x==0.0 && y<0.0
#            phi= 3*pi/2
#    elseif x==0.0 && y==0.0
#            phi= 0.0
#    elseif y == 0.0 && x < 0.0
#            phi = pi
#    elseif y == 0.0 && x > 0.0
#            phi = -pi
#    else
            phi = atan(y,x)
#    end

    return [R theta phi]
end





function SepVec(i,j,SP::SystemParameters) 
    ri = SP.BaseVecs[i,:]
    rj = SP.BaseVecs[j,:]
    
    return ri-rj
end

###################################
## SPherical Harmonic Y_l^m(theta,phi)
##################################
function SphY(l,m,theta,phi)

   l=Int(l)
   m=Int(m)

   if(m < 0 )
     return (-1.0)^(abs(m))*sf_legendre_sphPlm(l, abs(m),cos(theta))*(cos(abs(m)*phi)-im*sin(abs(m)*phi))
   else
     return sf_legendre_sphPlm(l, m,cos(theta))*(cos(m*phi)+im*sin(m*phi))
   end #endif

end
#############################################
## Derivatives with respect to theta and Phi
##############################################

function DthSphY(l,m,th,ph)

    if l!=m
        if (th == 0.0 || th == pi)
         return sqrt((l-m)*(l+m+1))SphY(l,m+1,th,ph)*(cos(ph)-im*sin(ph))
        else
      
         return m*SphY(l,m,th,ph)*cot(th)+sqrt((l-m)*(l+m+1))SphY(l,m+1,th,ph)*(cos(ph)-im*sin(ph))
      
       end
    
    else 
        return 0.0
    end
        
        
end
########### Dphi  ################
function DphiSphY(l,m,theta,phi)
    return im*m*SphY(l,m,theta,phi)
end
#################################################

 ###################################
 ## 2da Deriviada con respecto a theta de
 ## SPherical Harmonic Y_l^m(theta,phi)
 ##################################

function DDSphY(l,m,th,ph)

     if (th == 0.0 || th == pi)
      return 0.0
     else
      return m*(m*(cot(th))^2 - (csc(th))^2)*SphY(l,m,th,ph)
      +sqrt((l-m)*(l+m+1)*(2*m+1))*cot(th)*SphY(l,m+1,th,ph)*(cos(ph)-im*sin(ph))
      +sqrt((l-m)*(l-m-1)*(l+m+2)*(l+m+1))*SphY(l,m+2,th,ph)*(cos(2*ph)-im*sin(2*ph))
    end #if
end #function

###################################
## Depolarization factor of sphere
##################################
N0(l) =  l/(2*l+1)


###################################
## Defines Kroncker Delta
##################################
KD(l1,l2) = ifelse(l1 == l2 , 1 , 0 )

function Aelement(l1, l2, m1, m2, i,j,SP::SystemParameters)
    
    Rij    = CartToSphere( SepVec(i,j,SP::SystemParameters) )[1]
    theta  = CartToSphere( SepVec(i,j,SP::SystemParameters) )[2]
    phi    = CartToSphere( SepVec(i,j,SP::SystemParameters) )[3]
    
    SPH = (-1.0)^(m2)*conj(SphY(l1+l2,m1-m2,theta,phi))
    
    Num1   = (4*pi)^3*factorial( BigInt(l1+l2+m1-m2) )*factorial( BigInt(l1+l2-m1+m2) )
    Den1   = factorial(BigInt(l1+m1))*factorial( BigInt(l1-m1) )*factorial(BigInt(l2+m2))*factorial( BigInt(l2-m2) )
    Den2   = (2*l1+1)*(2*l2+1)*(2l1+2l2+1)
    
    if i==j
        return 0.0
    else
        return Complex{Float64}(sqrt(Num1/(Den1*Den2))* (SPH/Rij^(l1+l2+1) )   ) 
    end
    
end #endfunction


function Helement(Indx1::Int,Indx2::Int,lmax::Int,SP::SystemParameters)
    
    MatVars = SP.IndxF(Indx1,Indx2,lmax,SP)
    l1      = MatVars[1]
    l2      = MatVars[2]
    m1      = MatVars[3]
    m2      = MatVars[4]
    i       = MatVars[5]
    j       = MatVars[6]

    a = SP.Radius
    
    term1 = N0(l1)*KD(l1,l2)*KD(m1,m2)*KD(i,j)
    term2 = sqrt(l1*l2*a[i]^(2*l1+1)*a[j]^(2*l2+1))/(4*pi)
    term3 = Aelement(l1,l2,m1,m2,i,j,SP::SystemParameters)
    
    return  Complex{Float64}(term1 + (-1.0)^l2*term2*term3)

end #endfunction

function Hmat(lmax::Int,SP::SystemParameters) 
    
    MatDim = SP.IndxF(0,0,lmax,SP)
    
    HH     = Complex{Float64}[Helement(Indx_i,Indx_j,lmax,SP::SystemParameters) for Indx_i = 1:MatDim, Indx_j = 1:MatDim]
    
    return HH
end
    
##--------------------------------------------------------------------
##    Periodic Interaction Matrices
##--------------------------------------------------------------------



function Aelement_periodic(l1, l2, m1, m2, i,j,mm,nn,k0vec,SP::SystemParameters)
    ri     = SP.BaseVecs[i,:]
    rj     = SP.BaseVecs[j,:]
    t1     = SP.LattVecs[1,:]
    t2     = SP.LattVecs[2,:]

    rij_mn = (rj + mm*t1+nn*t2) - ri 
    
     Rij    = InteractionMatrix.CartToSphere( rij_mn )[1]
    theta  = InteractionMatrix.CartToSphere( rij_mn )[2]
    phi    = InteractionMatrix.CartToSphere(  rij_mn )[3]
     SPH = (-1.0)^(m2)*conj(InteractionMatrix.SphY(l1+l2,m1-m2,theta,phi))
    
    Num1   = (4*pi)^3*factorial( BigInt(l1+l2+m1-m2) )*factorial( BigInt(l1+l2-m1+m2) )
    Den1   = factorial(BigInt(l1+m1))*factorial( BigInt(l1-m1) )*factorial(BigInt(l2+m2))*factorial( BigInt(l2-m2) )
    Den2   = (2*l1+1)*(2*l2+1)*(2l1+2l2+1)
    
    if Rij == 0 #i==j && mm == 0 && nn == 0
        return 0.0
    else
        out = Complex{Float64}(sqrt(Num1/(Den1*Den2))* (SPH/Rij^(l1+l2+1) )   )*exp(-im*dot(k0vec,rij_mn)) 
        
        return out
    end
    
end #endfunction


KD2(i,j) = ifelse(i == j , 0 , 1 )
KD3(x,x0) = ifelse(x <= x0 + 0.04 , 1 , 0 )

function KD4(x,x0) 
    
    if x[3] != 0
        return 1
    else
        return KD3(norm(x[1:2]),x0) 
    end
    
end

function Helement_periodic(Indx1::Int,Indx2::Int,NCells,lmax::Int,k0vec,SP::SystemParameters)
    
    LVecs  = SP.LattVecs
    Rs     = SP.BaseVecs
    a      = SP.Radius    
    
    MatVars = SP.IndxF(Indx1,Indx2,lmax,SP)
    l1      = MatVars[1]
    l2      = MatVars[2]
    m1      = MatVars[3]
    m2      = MatVars[4]
    i       = MatVars[5]
    j       = MatVars[6]

    
    term1 = InteractionMatrix.N0(l1)*KD(l1,l2)*KD(m1,m2)*KD(i,j)
    term2 = sqrt(l1*l2*a[i]^(2*l1+1)*a[j]^(2*l2+1))/(4*pi)
    term3 = 0.0+im*0.0
    
        
    for mm = -NCells:NCells, nn = -NCells:NCells
        ri     = SP.BaseVecs[i,:]
        rj     = SP.BaseVecs[j,:]
        t1     = SP.LattVecs[1,:]
        t2     = SP.LattVecs[2,:]

        rij_mn = (rj + mm*t1+nn*t2) - ri 
        
        term3 += Aelement_periodic(l1,l2,m1,m2,i,j,mm,nn,k0vec,SP::SystemParameters)*KD3(norm(rij_mn),SP.NNSep)
    end
        

    
    out_r = real(Complex{Float64}(term1 + (-1.0)^l2*term2*term3))
    out_i = imag(Complex{Float64}(term1 + (-1.0)^l2*term2*term3))
    norm(out_r) <= 1e-8 ? out_r = 0.0 : out_r = out_r
    norm(out_i) <= 1e-8 ? out_i = 0.0 : out_i = out_i
  
    return out_r+im*out_i

end #endfunction


function InteractionSum_H(Indx1,Indx2,lmax,NCells,k0vec,NInteractions,SP::SystemParameters)
    
    MatVars = SP.IndxF(Indx1,Indx2,lmax,SP)
    Nint    = NInteractions
    l1      = MatVars[1]
    l2      = MatVars[2]
    m1      = MatVars[3]
    m2      = MatVars[4]
    i       = MatVars[5]
    j       = MatVars[6]
    MatDim  = SP.IndxF(0,0,lmax,SP)
    NNCSep  = SP.NNSep
    
    Dims =  0
    
    term3 = 0.0+im*0.0
    t1    = SP.LattVecs[1,:]
    t2    = SP.LattVecs[2,:]
    
    norm(t2) == 0 ? Dims = 1 : Dims = 2
    
   
    if Dims == 2 && SP.NP > 1
        

        for mm = -NCells:NCells, nn = -NCells:NCells
            ri     = SP.BaseVecs[i,:]
            rj     = SP.BaseVecs[j,:]

            rij_mn = (rj + mm*t1+nn*t2) - ri 
            #InPlaneCSep = norm( (rj[1:2] + mm*t1[1:2]+nn*t2[1:2]) - ri[1:2] ) 
            term3 += Aelement_periodic(l1,l2,m1,m2,i,j,mm,nn,k0vec,SP::SystemParameters)*KD4(rij_mn,Nint*NNCSep)
        end
        
    elseif Dims == 2 && SP.NP ==1
        

        for mm = -NCells:NCells, nn = -NCells:NCells
            ri     = SP.BaseVecs[i,:]

            rij_mn = ( mm*t1+nn*t2) - ri 
            InPlaneCSep = norm( mm*t1[1:2]+nn*t2[1:2] - ri[1:2] ) 

            term3 += Aelement_periodic(l1,l2,m1,m2,i,j,mm,nn,k0vec,SP::SystemParameters)*KD3(InPlaneCSep,NNCSep)
        end
        
    elseif Dims == 1 && SP.NP > 1
         
        for mm = -NCells:NCells
            
            nn = 0
            ri = SP.BaseVecs[i,:]
            rj = SP.BaseVecs[j,:]

            rij_mn = (rj + mm*t1+nn*t2) - ri 
            InPlaneCSep = norm( (rj[1:2] + mm*t1[1:2]+nn*t2[1:2]) - ri[1:2] )  
            term3 += Aelement_periodic(l1,l2,m1,m2,i,j,mm,nn,k0vec,SP::SystemParameters)*KD3(InPlaneCSep,NNCSep)
            
        end
        
        
    elseif Dims == 1 && SP.NP == 1
       
       for mm = -NCells:NCells
            
            nn = 0
            ri = SP.BaseVecs[i,:]
            

            rij_mn = (mm*t1+nn*t2)-ri 
             InPlaneCSep = norm( (mm*t1[1:2]+nn*t2[1:2]) - ri[1:2] )  
            term3 += Aelement_periodic(l1,l2,m1,m2,i,j,mm,nn,k0vec,SP::SystemParameters)*KD3(InPlaneCSep,NNCSep)
        end

    end
    
    return term3
    
end

function Helement_periodic_NearestNeighbors(Indx1::Int,Indx2::Int,NCells,lmax::Int,k0vec,Nint::Int,SP::SystemParameters)
    
    LVecs  = SP.LattVecs
    Rs     = SP.BaseVecs
    a      = SP.Radius        
    
    MatVars = SP.IndxF(Indx1,Indx2,lmax,SP)
    l1      = MatVars[1]
    l2      = MatVars[2]
    m1      = MatVars[3]
    m2      = MatVars[4]
    i       = MatVars[5]
    j       = MatVars[6]

    
    term1 = InteractionMatrix.N0(l1)*KD(l1,l2)*KD(m1,m2)*KD(i,j)
    term2 = sqrt(l1*l2*a[i]^(2*l1+1)*a[j]^(2*l2+1))/(4*pi)
    term3 = InteractionSum_H(Indx1,Indx2,lmax,NCells,k0vec,Nint,SP)
    
    out_r = real(Complex{Float64}(term1 + (-1.0)^l2*term2*term3))
    out_i = imag(Complex{Float64}(term1 + (-1.0)^l2*term2*term3))
    norm(out_r) <= 1e-9 ? out_r = 0.0 : out_r = out_r
    norm(out_i) <= 1e-9 ? out_i = 0.0 : out_i = out_i
  
    return out_r+im*out_i

end #endfunction


function Hmat_periodic(lmax::Int,NCells,k0vec,NInteractions,SP::SystemParameters) 
    
    MatDim = SP.IndxF(0,0,lmax,SP)
    Nint   = NInteractions
    #HH     = Complex[Helement_periodic(Indx_i,Indx_j,NCells,lmax::Int,k0vec,SP::SystemParameters) for Indx_i = 1:MatDim, Indx_j = 1:MatDim]
    HH     = Complex[Helement_periodic_NearestNeighbors(Indx_i,Indx_j,NCells,lmax::Int,k0vec,Nint,SP::SystemParameters) for Indx_i = 1:MatDim, Indx_j = 1:MatDim]
    
    return HH
end

# ############################################# 
# temp. removes interaction between layers
# ###########################################
function Hmat_periodic_temp(lmax::Int,NCells,k0vec,Inter,NInteractions,SP::SystemParameters) 
    
    MatDim = SP.IndxF(0,0,lmax,SP)
    Nint   = NInteractions
    #HH     = Complex[Helement_periodic(Indx_i,Indx_j,NCells,lmax::Int,k0vec,SP::SystemParameters) for Indx_i = 1:MatDim, Indx_j = 1:MatDim]
    HH     = Complex[Helement_periodic_NearestNeighbors_temp(Indx_i,Indx_j,NCells,lmax::Int,k0vec,Inter,Nint,SP::SystemParameters) for Indx_i = 1:MatDim, Indx_j = 1:MatDim]
    
    return HH
end

function particle_layer(i,particles_per_layer)
    
    if i <= particles_per_layer
        return 1.0
    else
        return 2.0
    end
    
end

function 

###########################
# Calculates a part of the interaction matrix.
## Only for bilayers with lmax = 1    
## Inter =  1 Interlayer (dif. layers)
## Inter = 2  Intralayer (same layuer)
## Only for monolayer with lmax > 1    
## Inter = 3  Dipolar   
##  Inter = 4 Multipolar
############################    
    Helement_periodic_NearestNeighbors_temp(Indx1::Int,Indx2::Int,NCells,lmax::Int,k0vec,Inter,Nint,SP::SystemParameters)
    LVecs  = SP.LattVecs
    Rs     = SP.BaseVecs
    a      = SP.Radius        
    
    MatVars = SP.IndxF(Indx1,Indx2,lmax,SP)
    l1      = MatVars[1]
    l2      = MatVars[2]
    m1      = MatVars[3]
    m2      = MatVars[4]
    i       = MatVars[5]
    j       = MatVars[6]
    NBase   = length(a)
    ppl     = Int(round(NBase/2) )
    
    layer_i     = particle_layer(i,ppl)
    layer_j     = particle_layer(j,ppl)
    
    term1 = InteractionMatrix.N0(l1)*KD(l1,l2)*KD(m1,m2)*KD(i,j)
    term2 = sqrt(l1*l2*a[i]^(2*l1+1)*a[j]^(2*l2+1))/(4*pi)
    term3 = InteractionSum_H(Indx1,Indx2,lmax,NCells,k0vec,Nint,SP)
    
    out_r = real(Complex{Float64}(term1 + (-1.0)^l2*term2*term3))
    out_i = imag(Complex{Float64}(term1 + (-1.0)^l2*term2*term3))
    norm(out_r) <= 1e-9 ? out_r = 0.0 : out_r = out_r
    norm(out_i) <= 1e-9 ? out_i = 0.0 : out_i = out_i

    # in-plane - out-plane (interlayer)
    if Inter == 1
        lmax != 1 ? error("for Inter = 1, lmax should be equal to 1")  :  i = i 
        #if abs(m1) != abs(m2) &&  3 < i+j <7
        if abs(m1) != abs(m2) &&  layer_i != layer_j
            return out_r+im*out_i   
        else
            return 0.0
        end
    # in-plane - in-plane (interlayer)
    elseif Inter == 2
        lmax != 1 ? error("for Inter = 2, lmax should be equal to 1")  :  i = i 

        if abs(m1) == abs(m2)  &&  layer_i != layer_j
            return out_r+im*out_i   
        else
            return 0.0
        end
    # total interlayer interaction
        
    elseif Inter == 0
        if   layer_i != layer_j
            return out_r+im*out_i   
        else
            return 0.0
        end
   
    #Intralayer (bilayers)
    elseif Inter == 3
        lmax != 1 ? error("for Inter = 2, lmax should be equal to 1")  :  i = i 

        if  layer_i == layer_j && i != j
            return out_r+im*out_i   
        else
            return 0.0
        end
    #dipole-quad and quad-quad (quadrupole)        
    elseif Inter == 4
        lmax == 1 ? error("for Inter = 4, lmax should be > 1")  :  i = i 

        if l1 != 1 || l2 != 1 
            return out_r+im*out_i   
        else
            return 0.0
        end
    #dipole-dipole (quadrupole approx.)        
        
     elseif Inter == 5
        lmax == 1 ? error("for Inter = 5, lmax should be > 1")  :  i = i 

         if l1 == 1 && l2 == 1 
            return out_r+im*out_i   
        else
            return 0.0
        end
        
    end
    
    
end #endfunction
###################################
# end temp
# ################################


## ------------------------------------------------ 
## Allows you to graph Hmat at varying kvec 
## for interaction between particles "i" and "j"
## -------------------------------------------------
function Hmat_PiPj(i,j,kmesh,lmax,SP::SystemParameters)
    
    MatDim   = SP.IndxF(0,0,lmax,SP)
    MatDimPP = Int(MatDim/SP.NP)
    indx_i1  = MatDimPP*(i-1) + 1
    indx_i2  = MatDimPP*i
    indx_j1  = MatDimPP*(j-1) + 1
    indx_j2  = MatDimPP*j
    NNmesh   = size(kmesh)[1]
    out   = []
    label = []
    count = 1
    for ii = 1 : NNmesh
        HH      = Complex[InteractionMatrix.Helement_periodic_NearestNeighbors(i,j,1,lmax,kmesh[ii,1:3],SP) for i = indx_i1:indx_i2, j = indx_j1:indx_j2]
        for i = 1:MatDimPP, j = 1 : MatDimPP
            j >= i ? out = push!(out,HH[i,j] ) : j = j
            
            if j >= i && ii == 1
                it = SP.IndxF(i,j,lmax,SP)
                #ts = string("(",it[1],it[2],it[3],it[4],")")
                ts = string("(",it[3],it[4],")")                
                count == 1 ? label = ts : label = hcat(label,ts)
                count +=1
            end
        end
        
    end
    
    return (Array( transpose(reshape(out,:,NNmesh) ) ), label ) 

end

end
