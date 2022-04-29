module SolveSystem

using  DielectricFunctions
using  InitializeSystem
import InteractionMatrix
using  LinearAlgebra
import PeriodicSystems

export Flm
export SolveDipoleAtLambda
export SolveCabsSpectrum
export ExtPlaneWavePotential
export SphToCart_Dipole,SolveCabs_2Dperiodic,SolveQs_Periodic
export ExtractSphericalDipole,SolveCabsSpectrum_AtKpoint
export SolveCabs_2Dperiodic_t,SolveCabsSpectrum_AtKpoint_t

function SolveQs_Periodic(lambda,lmax ,polarization ,k0vec ,Nint,SP::SystemParameters)
    
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    Lvecs  = SP.LattVecs
    CellArea = norm(Lvecs[1,:])^2
    
    MatDim = SP.IndxF(0,0,lmax,SP)   
#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat_periodic(lmax::Int,Nint,k0vec,Nint,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)
    U2  = transpose(conj.(U)) 
    
    function Gelement(i,j,lambda)
        
        ii = SP.IndxF(i,j,lmax,SP)[6]
        
        Gij=0.0
        for s=1:MatDim
            Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    E0Sph   =  SolveSystem.ExtPlaneWavePotential(lmax,polarization,SP).*Phase(k0vec,lmax,SP::SystemParameters)
       
    GG = [Gelement(i,j,lambda) for i = 1:MatDim, j = 1:MatDim]
        
    Qs     = GG*E0Sph
     
    for i = 1:length(Qs)
        ll = SP.IndxF(1,i,lmax,SP)[2]
        j = SP.IndxF(1,i,lmax,SP)[6]
        Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
     end
        
      
  return Qs
    
end

function Phase(kvec,lmax,SP::SystemParameters)
    re1 = []
    im1 = []
    MatDim = SP.IndxF(0,0,lmax,SP)
    Rs = SP.BaseVecs
    for i = 1: MatDim
        j  =  SP.IndxF(i,1,lmax,SP)[5]
        re1 = push!(re1, real( exp(im*dot(kvec,Rs[j,:])) ) )    
        im1 = push!(im1, imag( exp(im*dot(kvec,Rs[j,:])) ) )    

    end

    return re1.+im*im1
    
end

function Flm(m1,a)

  if m1==0
    Fl = (sqrt(a^3)/(4*pi))*sqrt(4pi/3)
  elseif m1==1
    Fl = (sqrt(a^3)/(4*pi))*sqrt(2pi/3)#*(-1+im)
  elseif m1==-1
    Fl = (sqrt(a^3)/(4*pi))*sqrt(2pi/3)#*(1+im)
  else
    error(" m must be -1,0,1")
  end

return Fl

end



#-----------------------------------------------------
# See: "Substrate effects on the optical properties of spheroidal nanoparticles": roman, noguez,barrera
# polarization is a string with values "x" "y" and "z".
# --------------------------------------------------

function ExtPlaneWavePotential(lmax,polarization,SP::SystemParameters)
    
    MatDim = SP.IndxF(0,0,lmax,SP)    
    Rs     = SP.BaseVecs
    a      = SP.Radius
    E0r = []
    
    if polarization     == "x"
        
        for i = 1: MatDim
            
            ll =  SP.IndxF(i,1,lmax,SP)[1]
            mm =  SP.IndxF(i,1,lmax,SP)[3]
            j  =  SP.IndxF(i,1,lmax,SP)[5]
            
            if ll == 1 && mm == 1
                E0r = push!(E0r,-real(Flm(1,a[j])) )
                
            elseif ll == 1 && mm == -1
                E0r = push!(E0r,real(Flm(-1,a[j])) )
                
            else
                E0r = push!(E0r,0)
            end
            
        end

        return 0.5*E0r
        
    elseif polarization == "y"
        
        for i = 1: MatDim
            
            ll =  SP.IndxF(i,1,lmax,SP)[1]
            mm =  SP.IndxF(i,1,lmax,SP)[3]
            j  =  SP.IndxF(i,1,lmax,SP)[5] 
            
            if ll == 1 && mm == 1
                E0r = push!(E0r,real(Flm(1,a[j])))
                
            elseif ll == 1 && mm == -1
                E0r = push!(E0r,real(Flm(1,a[j])) )
                
            else
                E0r = push!(E0r,0)
            end
            
        end


        return im*0.5*(E0r)
        
    elseif polarization == "z"
        
        for i = 1: MatDim
            
            ll =  SP.IndxF(i,1,lmax,SP)[1]
            mm =  SP.IndxF(i,1,lmax,SP)[3]
            j  =  SP.IndxF(i,1,lmax,SP)[5]
            
            if ll == 1 && mm == 0
                E0r = push!(E0r,Flm(0,a[j]))
            else
                E0r = push!(E0r,0)
            end
            
        end

        return E0r
    else
        error("Polarization must be x, y or z.\n CHECK: SolveSystem.ExtPlaneWave ")
        
    end
    

end

function ExtPlaneWaveCart(polarization,SP::SystemParameters)

    for i = 1:SP.NP
        
        if polarization == "x"
          return [1,0,0]
            
        elseif polarization == "y"
          return [0,1,0]
        
        elseif polarization == "z"
          return [0,0,1]

        end
        
    end
end

function ExtPlaneWaveCartPeriodic(polarization,kvec,SP::SystemParameters)
    Rs = SP.BaseVecs
    Eout_r = []
    Eout_i = []
    
    for i = 1:SP.NP
    

        
        if polarization == "x"
            
          Eout_r = append!(Eout_r,real.([1,0,0]*exp(im*dot(kvec,Rs[i,:])) ) )
          Eout_i = append!(Eout_i,imag.([1,0,0]*exp(im*dot(kvec,Rs[i,:])) ) )
           
        elseif polarization == "y"
          Eout_r = append!(Eout_r,real.([0,1,0]*exp(im*dot(kvec,Rs[i,:])) ) )
          Eout_i = append!(Eout_i,imag.([0,1,0]*exp(im*dot(kvec,Rs[i,:])) ) )
        
        elseif polarization == "z"
          Eout_r = append!(Eout_r,real.([0,0,1]*exp(im*dot(kvec,Rs[i,:])) ) )
          Eout_i = append!(Eout_i,imag.([0,0,1]*exp(im*dot(kvec,Rs[i,:])) ) )

        end
        
    end
    
    Eout_r = transpose( reshape(Eout_r,3,:) )
    Eout_i = transpose( reshape(Eout_i,3,:) )
   
    return Array(Eout_r.+im*Eout_i)

end



function SolveQsAtLambda(lambda,lmax,polarization,SP::SystemParameters)
    
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    MatDim = SP.IndxF(0,0,lmax,SP)

#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat(lmax::Int,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)

    function Gelement(i,j)
        ii = SP.IndxF(i,j,lmax,SP)[5]
        Gij=0.0
        
        for s=1:lmax
            Gij-=(U[i,s]*U[j,s])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction
    
    E0Sph   =  -ExtPlaneWavePotential(lmax,polarization,SP)
       
    GG   = [Gelement(i,j) for i = 1:MatDim, j = 1:MatDim]
    Qs   = GG*E0Sph 
    
    for i = 1:length(Qs)
        ll = SP.IndxF(1,i,lmax,SP)[2]
        j  = SP.IndxF(1,i,lmax,SP)[6]
        Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))s
    end
     
    return Qs

end#endfunction


#############################################
## Note: returns dipole moments in spherical coordinates
##      for a given lambda and a given "m" (external potentail)
########################################

function SolveDipoleAtLambda(lambda,lmax,polarization,SP::SystemParameters)
    
    Qs     = SolveQsAtLambda(lambda,lmax,polarization,SP::SystemParameters)     
    PsSph  = ExtractSphericalDipole(Qs,lmax,SP::SystemParameters)
    Ps     = SphToCart_Dipole(PsSph,lmax,SP::SystemParameters)
            
  return Ps

end#endfunction

function SolveCabsSpectrum(;lambda_list="",lmax ="",polarization = "",SP::SystemParameters ="")
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    MatDim = SP.IndxF(0,0,lmax,SP)   
#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat(lmax::Int,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)

    function Gelement(i,j,lambda)
        
        ii = SP.IndxF(i,j,lmax,SP)[6]
        
        Gij=0.0
        for s=1:MatDim
            Gij-=(U[i,s]*U[j,s])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    E0Sph   =  ExtPlaneWavePotential(lmax,polarization,SP)
    E0      =  ExtPlaneWaveCart(polarization,SP)
   
    cabs_out = []
    
    for lambdai in lambda_list
        ki       = 2pi*sqrt(Emed)/lambdai
        GG = [Gelement(i,j,lambdai) for i = 1:MatDim, j = 1:MatDim]
        
        Qs     = GG*E0Sph
        
        for i = 1:length(Qs)
            ll = SP.IndxF(1,i,lmax,SP)[2]
            j = SP.IndxF(1,i,lmax,SP)[6]
            Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
        end
        
        PsSph  = -ExtractSphericalDipole(Qs,lmax,SP)
        Ps     = SphToCart_Dipole(PsSph,lmax,SP::SystemParameters)
        cabs_i = 0.0
        
        for i = 1:SP.NP
            cabs_i   += 4pi*ki*imag(dot(Ps[i],E0) )

        end
        
        cabs_out = push!(cabs_out,cabs_i)

    end
    
  return cabs_out

end#endfunction

function ExtractSphericalDipole(Qs,lmax,SP::SystemParameters)
        
    Psr   = []
    Psi   = [] 
    MatDim = SP.IndxF(0,0,lmax,SP)   
   
    for i = 1: MatDim
        
        ll =  SP.IndxF(i,1,lmax,SP)[1]
        mm =  SP.IndxF(i,1,lmax,SP)[3]
        j  =  SP.IndxF(i,1,lmax,SP)[5]
            
                        
            if ll == 1 && mm == 1
                Psr = push!(Psr,real(Qs[i]) )
                Psi = push!(Psi,imag(Qs[i])  )
                
            elseif ll == 1 && mm == -1
                Psr = push!(Psr,real(Qs[i]) )
                Psi = push!(Psi,imag(Qs[i])  )
                
            elseif ll == 1 && mm== 0
                Psr = push!(Psr,real(Qs[i]) )
                Psi = push!(Psi,imag(Qs[i])  )
                
            end
            
            
        end
    
    return Psr+im*Psi

end
    
function SphToCart_Dipole(Ps,lmax,SP::SystemParameters)
        
    Psr   = 0
    Psi   = 0
    MatDim = SP.IndxF(0,0,lmax,SP)   

    Pout =[] 
       
    for i = 1: SP.NP
        i1 = (i-1)*3+1
        i2 = (i-1)*3+2
        i3 = (i-1)*3+3
        
        px = sqrt( (8pi)/3 ) *(Ps[i1]-Ps[i3])/2                        
        py = sqrt( (8pi)/3 ) *(Ps[i1]+Ps[i3])/(2*im)
        pz = sqrt(4pi/3)*Ps[i2]
        
        Pout = push!(Pout,[px py pz])
            
        
    end
        
        return Pout

end

##--------------------------------------------------------------------
##   Periodic Functions
##--------------------------------------------------------------------

function SolveQsAtLambda(lambda,lmax,polarization,NCells,k0vec,SP::SystemParameters)
    
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    MatDim = SP.IndxF(0,0,lmax,SP)

#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat_periodic(lmax::Int,NCells,k0vec,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)

    function Gelement(i,j)
        ii = SP.IndxF(i,j,lmax,SP)[5]
        Gij=0.0
        
        for s=1:lmax
            Gij-=(U[i,s]*U[j,s])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction
    
    E0Sph   =  -SolveSystem.ExtPlaneWavePotential(lmax,polarization,SP)
       
    GG   = [Gelement(i,j) for i = 1:MatDim, j = 1:MatDim]
    Qs   = GG*(E0Sph.*Phase(k0vec,lmax,SP) )
    
    for i = 1:length(Qs)
        ll = SP.IndxF(1,i,lmax,SP)[2]
        j = SP.IndxF(1,i,lmax,SP)[6]
        Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
    end
     
    return Qs

end#endfunction



function SolveCabsSpectrum(;lambda_list="",lmax ="",polarization = "",k0vec = "",SP::SystemParameters ="")
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    Lvecs  = SP.LattVecs
    CellArea = norm(cross(Lvecs[1,:],Lvecs[2,:]))
    
    MatDim = SP.IndxF(0,0,lmax,SP)   
#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat_periodic(lmax::Int,1,k0vec,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)

    function Gelement(i,j,lambda)
        
        ii = SP.IndxF(i,j,lmax,SP)[6]
        
        Gij=0.0
        for s=1:MatDim
            Gij-=(U[i,s]*U[j,s])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    E0Sph   =  SolveSystem.ExtPlaneWavePotential(lmax,polarization,SP)
    E0      =  SolveSystem.ExtPlaneWaveCart(polarization,SP)
   
    cabs_out = []
    
    for lambdai in lambda_list
        ki       = 2pi*sqrt(Emed)/lambdai
        GG = [Gelement(i,j,lambdai) for i = 1:MatDim, j = 1:MatDim]
        
        Qs     = GG*E0Sph
        
        for i = 1:length(Qs)
            ll = SP.IndxF(1,i,lmax,SP)[2]
            j = SP.IndxF(1,i,lmax,SP)[6]
            Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
        end
     
        PsSph  = -SolveSystem.ExtractSphericalDipole(Qs,lmax,SP)
        Ps     = SolveSystem.SphToCart_Dipole(PsSph,lmax,SP::SystemParameters)
        cabs_i = 0.0
        
        for i = 1:SP.NP
            cabs_i   += 4pi*ki*imag(dot(Ps[i],E0) )/CellArea
        end
        
        cabs_out = push!(cabs_out,cabs_i)

    end
    
  return cabs_out

end#endfunction



function SolveCabsSpectrum_AtKpoint(lambda_list,lmax ,polarization ,k0vec ,NInteractions,SP::SystemParameters)
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    Lvecs  = SP.LattVecs
    CellArea = norm(cross(Lvecs[1,:],Lvecs[2,:]) )
    
    MatDim = SP.IndxF(0,0,lmax,SP)   
#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat_periodic(lmax::Int,1,k0vec,NInteractions,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)
    U2  = transpose(conj.(U)) 
    
    function Gelement(i,j,lambda)
        
        ii = SP.IndxF(i,j,lmax,SP)[6]
        
        Gij=0.0
        for s=1:MatDim
            #Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,wp1,gamma1,Emed,a[ii])-ns[s])
            Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    E0Sph =SolveSystem.ExtPlaneWavePotential(lmax,polarization,SP).*Phase(k0vec,lmax,SP::SystemParameters)
    E0      =  ExtPlaneWaveCartPeriodic(polarization,k0vec,SP::SystemParameters)
   
    cabs_out = []
    
    for lambdai in lambda_list
        ki = 2pi*sqrt(Emed)/lambdai
        GG = [Gelement(i,j,lambdai) for i = 1:MatDim, j = 1:MatDim]
        
        Qs     = GG*E0Sph
        
        for i = 1:length(Qs)
            ll    = SP.IndxF(1,i,lmax,SP)[1]
            j     = SP.IndxF(1,i,lmax,SP)[5]
            Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
        end
     
        PsSph  = -SolveSystem.ExtractSphericalDipole(Qs,lmax,SP)
  
        Ps     = SolveSystem.SphToCart_Dipole(PsSph,lmax,SP::SystemParameters)
        cabs_i = 0.0
        
        for i = 1:SP.NP
            cabs_i   += 4pi*ki*imag(dot(Ps[i],E0[i,:]) )/CellArea
        end
        
        cabs_out = push!(cabs_out,cabs_i)

    end
    
  return cabs_out

end#endfunction



function SolveCabs_2Dperiodic(eV_list, kmesh, lmax, polarization,Nint, SP::SystemParameters)
    
    dataOut = []
    Nkmesh  = size(kmesh)[1]
    lambda_list = lte.(eV_list)
    
    for i = 1:Nkmesh
        data_ki = SolveCabsSpectrum_AtKpoint(lambda_list,lmax ,polarization , kmesh[i,1:3] ,Nint,SP) 
        dataOut = append!(dataOut,data_ki)

    end
    
    return Array{Float64}(reshape( dataOut,  length(lambda_list),: ) )
    
end




function SolveCabsSpectrum_AtKpoint_ll(lambda_list,lmax ,polarization ,k0vec ,NInteractions,SP::SystemParameters)
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    Lvecs  = SP.LattVecs
    CellArea = norm(cross(Lvecs[1,:],Lvecs[2,:]) )
    
    MatDim = SP.IndxF(0,0,lmax,SP)   
#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat_periodic(lmax::Int,1,k0vec,NInteractions,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)
    U2  = transpose(conj.(U)) 
    
    function Gelement(i,j,lambda)
        
        ii = SP.IndxF(i,j,lmax,SP)[6]
        
        Gij=0.0
        for s=1:MatDim
            #Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,wp1,gamma1,Emed,a[ii])-ns[s])
            Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    E0Sph =SolveSystem.ExtPlaneWavePotential(lmax,polarization,SP).*Phase(k0vec,lmax,SP::SystemParameters)
    E0      =  ExtPlaneWaveCartPeriodic(polarization,k0vec,SP::SystemParameters)
   
    cabs_out = []
    
    for lambdai in lambda_list
        ki = 2pi*sqrt(Emed)/lambdai
        
        if norm(k0vec) > ki
            cabs_out = push!(cabs_out,0.0)
        else

            GG = [Gelement(i,j,lambdai) for i = 1:MatDim, j = 1:MatDim]
        
            Qs     = GG*E0Sph
        
            for i = 1:length(Qs)
                ll    = SP.IndxF(1,i,lmax,SP)[1]
                j     = SP.IndxF(1,i,lmax,SP)[5]
                Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
            end
     
            PsSph  = -SolveSystem.ExtractSphericalDipole(Qs,lmax,SP)
  
            Ps     = SolveSystem.SphToCart_Dipole(PsSph,lmax,SP::SystemParameters)
            cabs_i = 0.0
        
            for i = 1:SP.NP
                cabs_i   += 4pi*ki*imag(dot(Ps[i],E0[i,:]) )/CellArea
            end
        
            cabs_out = push!(cabs_out,cabs_i)
        end

    end
    
  return cabs_out

end#endfunction

function SolveCabs_2Dperiodic_ll(eV_list, kmesh, lmax, polarization,Nint, SP::SystemParameters)
    
    dataOut = []
    Nkmesh  = size(kmesh)[1]
    lambda_list = lte.(eV_list)
    
    for i = 1:Nkmesh
        data_ki = SolveCabsSpectrum_AtKpoint_ll(lambda_list,lmax ,polarization , kmesh[i,1:3] ,Nint,SP) 
        dataOut = append!(dataOut,data_ki)

    end
    
    return Array{Float64}(reshape( dataOut,  length(lambda_list),: ) )
    
end

## ##########################3
## temp stuff erase after
################################
function SolveCabsSpectrum_AtKpoint_t(lambda_list,lmax ,polarization ,k0vec ,SP::SystemParameters)
    Emed   = SP.Emed
    a      = SP.Radius
    Specie = SP.Specie
    Lvecs  = SP.LattVecs
    CellArea = SP.NP*pi*a[1]*a[1]#norm(cross(Lvecs[1,:],Lvecs[2,:]) )
    
    MatDim = SP.IndxF(0,0,lmax,SP)   
#Geometric function eval. at hw
    HH  = InteractionMatrix.Hmat_periodic(lmax::Int,1,k0vec,SP::SystemParameters)
    ns  = eigvals(HH)
    U   = eigvecs(HH)
    U2  = transpose(conj.(U)) 
    
    function Gelement(i,j,lambda)
        
        ii = SP.IndxF(i,j,lmax,SP)[6]
        
        Gij=0.0
        for s=1:MatDim
            #Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,wp1,gamma1,Emed,a[ii])-ns[s])
            Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    E0Sph =SolveSystem.ExtPlaneWavePotential(lmax,polarization,SP)#.*Phase(k0vec,lmax,SP::SystemParameters)
    E0      =  ExtPlaneWaveCartPeriodic(polarization,k0vec,SP::SystemParameters)
   
    cabs_out = []
    
    for lambdai in lambda_list
        ki = 2pi*sqrt(Emed)/lambdai
        GG = [Gelement(i,j,lambdai) for i = 1:MatDim, j = 1:MatDim]
        
        Qs     = GG*E0Sph
        
        for i = 1:length(Qs)
            ll    = SP.IndxF(1,i,lmax,SP)[1]
            j     = SP.IndxF(1,i,lmax,SP)[5]
            Qs[i] = Qs[i]*sqrt(ll*a[j]^(2ll+1))
        end
     
        PsSph  = -SolveSystem.ExtractSphericalDipole(Qs,lmax,SP)
        Ps     = SolveSystem.SphToCart_Dipole(PsSph,lmax,SP::SystemParameters)
        cabs_i = 0.0
        
        for i = 1:SP.NP
            cabs_i   += 4pi*imag(dot(Ps[i],E0[i,:]) )/CellArea

        end
        
        cabs_out = push!(cabs_out,cabs_i)

    end
    
  return cabs_out

end#endfunction



function SolveCabs_2Dperiodic_t(eV_list, kmesh, lmax, polarization, SP::SystemParameters)
    
    dataOut = []
    Nkmesh  = size(kmesh)[1]
    lambda_list = lte.(eV_list)
    
    for i = 1:Nkmesh
        data_ki = SolveCabsSpectrum_AtKpoint_t(lambda_list,lmax ,polarization , kmesh[i,1:3] ,SP) 
        dataOut = append!(dataOut,data_ki)

    end
    
    return Array{Float64}(reshape( dataOut,  length(lambda_list),: ) )
    
end


####
## end temp
#####
end