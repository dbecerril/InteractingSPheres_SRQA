module PeriodicSystems

using LinearAlgebra
using DelimitedFiles
using Distributed

using DielectricFunctions
using InitializeSystem
import InteractionMatrix
import ElectricFieldMod
import SolveSystem
import BandUnfolding
bu      = BandUnfolding


export DispersionRelation,OrderedDispersionRelation,Gmat_AtKpoint_lambda
export MakeKspaceMesh,ReciprocalVectors,EigenAtKpoint
export EigenSortFromFiles,OrderBandsFromFiles,WriteBandsToFile
export ConvergeVaryingLmax,omegaf2,SolveCabs_JuliaPlot,SolveCabsSpectrum_AtKpoint
export DensityOfStates,DensityOfStatesList,omegaf3,IIW_kpath,IIW_kpath_quad
##--------------------------------------------------------------------
##    Dispersion Relation Functions
##--------------------------------------------------------------------
function omegaf(ns,omegap,gamma,SP::SystemParameters)
    Emed  = SP.Emed
    
    omegap  = GetDrudeParameters(SP)[1]
    gamma   = GetDrudeParameters(SP)[2]
    Einf    = GetDrudeParameters(SP)[3]
    
    out   = zeros(length(ns)) + im*zeros(length(ns))
  
    
    for ii = 1:length(ns)
        ni =ns[ii]
        #term1   = im*gamma*(Emed*(1-ni)+ni)
        #term2   = sqrt(-(Emed*(ni-1)-ni)*(4*ni*omegap^2 - Emed*gamma^2 + (Emed-1)*ni*gamma^2 ))
        #term3   = 2*Emed*(ni-1)-2*ni
        #out[ii] = (term1-term2)/term3
        A_ii    = 1/(ni*( Einf -Emed)+Emed)
        GGamma  = 0.5*gamma
        out[ii] = sqrt(A_ii*omegap^2*ni - GGamma )-im*GGamma
    end
    
    return real.(out)
    
end

function omegaf2(ns,SP::SystemParameters)
    Emed  = SP.Emed
    
    omegap  = GetDrudeParameters(SP)[1]
    gamma   = GetDrudeParameters(SP)[2]

    ni =ns
    term1   = im*gamma*(Emed*(1-ni)+ni)
    term2   = sqrt(-(Emed*(ni-1)-ni)*(4*ni*omegap^2 - Emed*gamma^2 + (Emed-1)*ni*gamma^2 ))
    term3   = 2*Emed*(ni-1)-2*ni
    out = (term1-term2)/term3
    
    return real.(out)
    
end


function omegaf3(ns,SP::SystemParameters)

    Emed  = SP.Emed

    omegap  = GetDrudeParameters(SP)[1]
    gamma   = GetDrudeParameters(SP)[2]
    Einf    = GetDrudeParameters(SP)[3]


    A_ii    = 1/(ns*( Einf -Emed)+Emed)
    GGamma  = 0.5*gamma
    out     = sqrt(A_ii*omegap^2*ns - GGamma^2 )-im*GGamma

    return real.(out)

end


function DispersionRelation_atKparallel(kparallel,omegap,gamma,lmax,Nint,SP::SystemParameters)
    
    HH       = InteractionMatrix.Hmat_periodic(lmax,1,kparallel,Nint,SP::SystemParameters)
    ns       = eigvals(HH)
    
    omega_ns = omegaf(ns,omegap,gamma,SP::SystemParameters)
    
    return sort(omega_ns)
end

function DispersionRelation(kmesh,lmax,Nint,SP::SystemParameters)
    omegap = GetDrudeParameters(SP)[1]
    gamma  = GetDrudeParameters(SP)[2]
    out    = [] 
    MatDim = SP.IndxF(0,0,lmax,SP::SystemParameters)
    
    for i = 1:length(kmesh[:,1])
        temp = DispersionRelation_atKparallel(kmesh[i,1:3],omegap,gamma,lmax,Nint,SP::SystemParameters)
        out  = append!(out,temp)
    end
    
    return Array( transpose(reshape(out,MatDim,:))   )
    
end

function eigensort(Eigen1,Eigen2,IndxOld)
    
    eigvals1 = Eigen1.values[IndxOld]    
    eigvecs1 = Eigen1.vectors[:,IndxOld]
    
    eigvals2 = Eigen2.values    
    eigvecs2 = Eigen2.vectors

    Dim     = length(eigvals1)
    IndxOut = vec(collect(1:Dim)) 
    
    for i = 1: Dim
        
        Perm  = Float64[abs(dot(eigvecs1[:,i],eigvecs2[:,j])) for j = 1:Dim]
        tempi = findmax(Perm)
        
        IndxOut[i] = tempi[2] 
        
        
    end
    
    return IndxOut
    
end

function OrderedDispersionRelation(kmesh,lmax,SP::SystemParameters)
    
    omegap  = GetDrudeParameters(SP::SystemParameters)[1]
    gamma   = GetDrudeParameters(SP::SystemParameters)[2]
    MatDim  = SP.IndxF(0,0,lmax,SP)
    out     = zeros(length(kmesh[:,1]),MatDim)
    IndxOld = collect(1:MatDim)

    HH1      = InteractionMatrix.Hmat_periodic(lmax,1,kmesh[1,:],SP::SystemParameters)
    Eigen1   = eigen(HH1)
    omega_ns = omegaf(Eigen1.values,omegap,gamma,SP::SystemParameters)

    out[1,:]   = transpose(omega_ns)
    
    for i = 2:length(kmesh[:,1])
        
        HH2    = InteractionMatrix.Hmat_periodic(lmax,1,kmesh[i,:],SP::SystemParameters)
        Eigen2 = eigen(HH2)
        
        IndxNew  = eigensort(Eigen1,Eigen2,IndxOld)
        omega_ns = omegaf(Eigen2.values[IndxNew],omegap,gamma,SP::SystemParameters)
        out[i,:] = transpose(omega_ns)
        
        IndxOld = IndxNew
        Eigen1  = Eigen2
    
    end
    
    return out
    
end

function MakeKspaceMesh(t1,t2,NSteps)
        
    g1 = ReciprocalVectors(t1,t2)[1]
    g2 = ReciprocalVectors(t1,t2)[2]
    
    indx1 = range(-1,stop = 1, length = NSteps )
        
    kx = []
    ky = []
    kz = []    
    
    for mm in indx1, nn in indx1
        
        ki = mm*g1+nn*g2
        kx = push!(kx,ki[1])    
        ky = push!(ky,ki[2])
        kz = push!(kz,ki[3])   
            
    end

    return hcat(kx,ky,kz)

end



function ReciprocalVectors(t1,t2)

    t3 = [0,0,1]
    g1 = (2pi/dot(cross(t2,t3),t1))*cross(t2,t3);
    g2 = (2pi/dot(cross(t3,t1),t2))*cross(t3,t1);
    
    return (g1,g2)
end



function FormatFileNameFromKpoint(x ::Int)

    MaxNoDigits   = 7
    NoDigits      = length(digits(x,base = 10))
    MissingDigits = MaxNoDigits - NoDigits

    if NoDigits < MaxNoDigits

        FormattedNo = string("/",0)

        for i = 1:MissingDigits - 1

            FormattedNo  = string(FormattedNo,0)
        end

        FormattedNo = string(FormattedNo,x,".txt")

        return FormattedNo

    else

        return string("/",x,".txt")
    end

end


## For a system SP::SystemParameters,lmax and a kpoint[kx,ky,kz,::int]
## finds eigenvals and eigenvecs and writes to a file in DirLocation
## Files are written ass column1 (Eigenvalues) Columns 2:MatDim (eigenvector)

function EigenAtKpoint(kpoint,lmax,SP::SystemParameters)
    DirLocation = string(pwd(),"/Temp",SP.SystemName)

    omegap  = GetDrudeParameters(SP::SystemParameters)[1]
    gamma   = GetDrudeParameters(SP::SystemParameters)[2]
    MatDim  = SP.IndxF(0,0,lmax,SP)
    IndxOld = collect(1:MatDim)

    HH1      = InteractionMatrix.Hmat_periodic(lmax,1,kpoint[1:3],SP::SystemParameters)
    Eigen1   = eigen(HH1)
    omega_ns = real.(Eigen1.values)
    
    FileName = FormatFileNameFromKpoint( Int(kpoint[4]) )
    io       = open(string(DirLocation,FileName),"w")
    out      = hcat(omega_ns,transpose(real.(Eigen1.vectors) ) )
    out2     = hcat(omega_ns,transpose(imag.(Eigen1.vectors) ) )

    out      = vcat(out,out2)

    writedlm(io,out )
    
    close(io)
    
end

function ConvergeVaryingLmax(lmax0 ::Int,lmaxf::Int,tolerance, kpoint::Array, SP::SystemParameters)
        
    converged = false
    RunNo     = 1
    lmax_i    = 0
    criteria  = 0
    while converged == false && lmax_i < lmaxf
     
        lmax_i    = lmax0 + (RunNo-1)*2

        HH1      = InteractionMatrix.Hmat_periodic(lmax_i,1,kpoint[1:3],SP::SystemParameters)
        Eigen1   = sort(real.(eigvals(HH1)) )
        HH2      = InteractionMatrix.Hmat_periodic(lmax_i+2,1,kpoint[1:3],SP::SystemParameters)
        Eigen2   = sort(real.(eigvals(HH2)) )

        criteria = abs.(Eigen1-Eigen2[1:length(Eigen1)])./abs.(Eigen1)
        println("testing lmax $(lmax_i), MatDim = $(size(HH2))\n Criteria at $(maximum(criteria) )")
        
        maximum(criteria) < tolerance ? converged = true : converged = false
        RunNo += 1 
    end


    if converged == true
        println("System converged at lmax = $(lmax_i) with $(criteria) error ")
    else
        println("System not converged using lmax in the range : $(lmax0) - $(lmaxxf)")
    end
end


function EigenSortFromFiles(eigvals1,eigvecs1,eigvals2,eigvecs2,IndxOld)

    eigvals1 = eigvals1[IndxOld]
    eigvecs1 = eigvecs1[:,IndxOld]

    Dim     = length(eigvals1)
    IndxOut = vec(collect(1:Dim))

    for i = 1: Dim

        Perm  = Float64[abs(dot(eigvecs1[:,i],eigvecs2[:,j])) for j = 1:Dim]
        tempi = findmax(Perm)

        IndxOut[i] = tempi[2]


    end

    return IndxOut

end

function ExtractEigenFromFile(FileLocation::String, lmax::Int, SP::SystemParameters)
    io1   = open(FileLocation)
    data1 = readdlm(io1)
    close(io1)

    MatDim   = SP.IndxF(0,0,lmax,SP)
    EigVals1 = real.(data1[1:MatDim,1])
    REigVecs1 = transpose(data1[1:MatDim,2:MatDim+1])  
    IEigVecs1 = transpose(data1[MatDim+1:2MatDim,2:MatDim+1])

    return (EigVals1,REigVecs1+im*IEigVecs1)
end

function OrderBandsFromFiles(SP::SystemParameters,Nkpoints,lmax)
    DirLocation = string(pwd(),"/Temp",SP.SystemName)    
    
    omegap  = GetDrudeParameters(SP::SystemParameters)[1]
    gamma   = GetDrudeParameters(SP::SystemParameters)[2]
   
    FileList  = readdir(DirLocation)
    file1     = string(DirLocation,"/",FileList[1])
 
    Eigen1   = ExtractEigenFromFile(file1,lmax,SP::SystemParameters)
    EigVals1 = Eigen1[1]
    EigVecs1 = Eigen1[2]

    MatDim  = SP.IndxF(0,0,lmax,SP)
    IndxOld = collect(1:MatDim)
    out     = zeros(Nkpoints,MatDim)
    
    omega_ns = omegaf(EigVals1[IndxOld],omegap,gamma,SP::SystemParameters)
    out[1,:] = transpose(omega_ns)
    
    for i = 2 : length(FileList) 
   
   
        file2 = string(DirLocation,"/",FileList[i])
        
        Eigen2   = ExtractEigenFromFile(file2,lmax,SP::SystemParameters)
        EigVals2 = Eigen2[1]
        EigVecs2 = Eigen2[2]
        
        IndxNew  = EigenSortFromFiles(EigVals1,EigVecs1,EigVals2,EigVecs2,IndxOld)
        omega_ns = omegaf(EigVals2[IndxNew],omegap,gamma,SP::SystemParameters)
        out[i,:] = transpose(omega_ns)

        IndxOld   = IndxNew
        EigVals1  = EigVals2
        EigVecs1  = EigVecs2

    end

    return out

    
end

function WriteBandsToFile(FileName,klist,bands,SP::SystemParameters) 
    file1    = string(pwd(),"/OutputFiles/",FileName)
    io       = open(file1,"w")
    writedlm(io,hcat(klist,bands))
    close(io)
    
end

wp1 = wpAg
gamma1 = 0.01*gammaAg

function Gmat_AtKpoint_lambda(i,j,lambda,k0vec ,lmax ,SP::SystemParameters)
    Emed   = SP.Emed
    a      = SP.Radius
    Rs     = SP.BaseVecs
    Specie = SP.Specie
    Lvecs  = SP.LattVecs
    CellArea = norm(Lvecs[1,:])^2
    
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
            Cs  = U[i,s]*U2[s,j]
            den = (uspec(lambda,Specie,Emed,a[ii])-ns[s])
            println("Cs = $(Cs) \n")
            println("den = $(den) \n")
            Gij-=(U[i,s]*U2[s,j])/(uspec(lambda,Specie,Emed,a[ii])-ns[s])
        end #endfor

        return Gij
    end #endfunction

    #GG = [Gelement(i,j,lambda) for i = 1:MatDim, j = 1:MatDim]
    
    return Gelement(i,j,lambda)

end#endfunction

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


#function SolveCabs_JuliaPlot(lambda_list, k_list, lmax, polarization, SP::SystemParameters)
#    
#    dataOut = []
#    NK  = length(k_list[:,1])
    
#    for i = 1: NK
#        data_ki = SolveCabsSpectrum_AtKpoint(lambda_list,lmax ,polarization , k_list[i,:] ,SP) 
#        dataOut = append!(dataOut,data_ki)

#    end
    
#    return reshape( dataOut, :, length(lambda_list) )
    
#end


function Lorentzian(x,x0,Gamma)
    term1 = 1/pi
    term2 = 0.5*Gamma
    term3 = (x-x0)^2+(0.5*Gamma)^2
    return term1*(term2/term3)
end

function DensityOfStates(x,E0List,CountList,Gamma)
    
    sum1 = 0
    
    for i = 1:length(E0List)
        sum1 += CountList[i]*Lorentzian(x,E0List[i],Gamma)
    end
    
    return sum1
end

function DensityOfStatesList(EigenEnergy,DeltaEnergy)
    
    Emax   = maximum(EigenEnergy)
    Emin   = minimum(EigenEnergy)
    Elist  = range(Emin,stop = Emax,step = DeltaEnergy)
    Evec   = sort(vec(EigenEnergy))
    EEN    = length(Elist)
    EOut   = []
    Ecount = [] 
    
    jstart = 1
    
    
    for i = 2:EEN
        count_i  = 1
       
        for j = jstart : length(Evec)
            
            if Elist[i-1]<= Evec[j] < Elist[i] 
                count_i += 1
            else
                EOuti   = Elist[i]+DeltaEnergy
                EOut    = push!(EOut,EOuti)
                Ecount  = push!(Ecount,count_i)
                jstart  = j
                break
            end
            
        end
        
    end
    
    return hcat(EOut,Ecount)

    
end

function cwns_kpath(Kpath,lmax,SP::SystemParameters)
    
    NN     = size(Kpath)[1]
    Kpath  = hcat(Kpath ,collect(1:NN) )
    dataCS = 0
    for i = 1: NN
        datai  = BandUnfolding.cw_ns(Kpath[i,:],lmax,system);
        i == 1 ? dataCS = datai : dataCS = vcat(dataCS,datai)
    end

    return dataCS
end

function InterlayerInteractionWeights2(kpoint,Inter,lmax,Nint,BL::SystemParameters)
    MatDim          = BL.IndxF(0,0,lmax,BL)
    DimsPerParticle = Int(MatDim/BL.NP)
    dgap = norm(BL.LattVecs[1,:])/sqrt(3)-2*BL.Radius[1]

    
    HH0      = InteractionMatrix.Hmat_periodic(lmax,1,kpoint[1:3],Nint,BL::SystemParameters)
    IIH      = InteractionMatrix.Hmat_periodic_temp(lmax,1,kpoint[1:3],Inter,Nint,BL::SystemParameters)
    EigenBL   = eigen(HH0)
    UU      = EigenBL.vectors
    ns      = real.(EigenBL.values)
        
    IIW    = []
        
    # --------------------------------------------- 
    # --- Runs over the eigenstates at point k
    # ---------------------------------------------
    for i = 1: MatDim
        wi  = dot( UU[:,i], (IIH)*UU[:,i] )/ns[i]  
        IIW = append!(IIW, wi)           
    end
    
    kindx  = kpoint[4]
    ks     = fill(norm(kindx),MatDim)
    ff(x)  = bu.omegaf3(x,BL)

    omegas = ff.(ns)
    return hcat(ks,real.(omegas),real.(IIW) )

end

# This one measures dip and dip-quad/quad-quad contributions
function InterlayerInteractionWeights3(kpoint,Nint,BL::SystemParameters)
    lmax            = 2
    MatDim          = BL.IndxF(0,0,lmax,BL)
    DimsPerParticle = Int(MatDim/BL.NP)
    dgap = norm(BL.LattVecs[1,:])/sqrt(3)-2*BL.Radius[1]

    
    HH0        = InteractionMatrix.Hmat_periodic(lmax,1,kpoint[1:3],Nint,BL::SystemParameters)
    Hdip       = InteractionMatrix.Hmat_periodic_temp(lmax,1,kpoint[1:3],5,Nint,BL::SystemParameters)
    Hdip_quad  = InteractionMatrix.Hmat_periodic_temp(lmax,1,kpoint[1:3],4,Nint,BL::SystemParameters)

    EigenBL   = eigen(HH0)
    UU      = EigenBL.vectors
    ns      = real.(EigenBL.values)
        
    IIW    = []
        
    # --------------------------------------------- 
    # --- Runs over the eigenstates at point k
    # ---------------------------------------------
    for i = 1: MatDim
        wi  = dot( UU[:,i], (Hdip_quad-Hdip)*UU[:,i] )/ns[i]  
        IIW = append!(IIW, wi)    
    end
    
    kindx  = kpoint[4]
    ks     = fill(norm(kindx),MatDim)
    ff(x)  = bu.omegaf3(x,BL)

    omegas = ff.(ns)
    return hcat(ks,real.(omegas),real.(IIW) )

end
# Calculates IP-IP/OP-OP and IP-OP
# Inter = 1  : in-plane - out-plane (interlayer)
# Inter == 2 : in-plane - in-plane (interlayer)
# Inter == 3 :  Intralayer (bilayers)
# Inter == 4 :  dipole-quad and quad-quad (quadrupole)     
# Inter == 5 : dipole-dipole (quadrupole approx.)  

function IIW_kpath(Kpath,Inter,lmax,Nint,SP::SystemParameters)
    
    NN     = size(Kpath)[1]
    Kpath  = hcat(Kpath ,collect(1:NN) )
    dataCS = 0
    for i = 1: NN
        datai  = InterlayerInteractionWeights2(Kpath[i,:],Inter,lmax,Nint,SP::SystemParameters)
        i == 1 ? dataCS = datai : dataCS = vcat(dataCS,datai)
    end

    return dataCS
end

function IIW_kpath_quad(Kpath,Inter,lmax,Nint,SP::SystemParameters)
    
    NN     = size(Kpath)[1]
    Kpath  = hcat(Kpath ,collect(1:NN) )
    dataCS = 0
    for i = 1: NN
        datai  = InterlayerInteractionWeights3(Kpath[i,:],Nint,SP::SystemParameters)
        i == 1 ? dataCS = datai : dataCS = vcat(dataCS,datai)
    end

    return dataCS
end




end