module ElectricFieldMod

using GSL
using LinearAlgebra
using GeometricalPredicates

using  DielectricFunctions
using  InitializeSystem
import InteractionMatrix
import SolveSystem


export PlotEF,PlotEF_periodic
export CartToSphere, ECartesian_periodic
export LSVFtoGCVF_periodic,PlotEF_2Dperiodic,MakeMesh_EF2D

function LSVFtoGCVF(i,Ar,Atheta,Aphi,x,y,z,SP::SystemParameters)
    
    Rs = SP.BaseVecs[i,:]
    
    x0 = Rs[1]
    y0 = Rs[2]
    z0 = Rs[3]
    
    z    = z-z0
    x    = x-x0
    y    = y-y0
    
    thet = CartToSphere([x y z])[2]
    phi  = CartToSphere([x y z])[3]


    Ax =Ar*sin(thet)*cos(phi) + Atheta*cos(thet)*cos(phi) - Aphi*sin(phi)
    Ay =Ar*sin(thet)*sin(phi) + Atheta*cos(thet)*sin(phi) + Aphi*cos(phi)
    Az =Ar*cos(thet)          - Atheta*sin(thet)

    Ax = Complex{Float64}(Ax)
    Ay = Complex{Float64}(Ay)
    Az = Complex{Float64}(Az)
   
    return [Ax, Ay, Az]
end

########################################################
### Goes from particle (local) Cartesian system to 
### Particle spherical system
### input :  x,y,z_i          (local)
### output:  R, thet, phi   (local)
###############################################
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
    #if x==0.0 && y==0.0
    #    if z > 0.0
    #        theta = 0.0
    #    elseif z < 0.0
    #        theta = pi
    #    end
    #elseif z>0
    #    theta = acos(z/R)
    #else 
    #    theta = acos(-z/R)
    #end
    theta = acos(z/R)
###### Defining phi

    phi = atan(y,x)
  

    return [R theta phi]
end

function InParticle(rvec,SP ::SystemParameters)

    Rs = SP.BaseVecs
    a  = SP.Radius

    out = false

    for i = 1 : length(Rs[:,1])
        norm(rvec-Rs[i,:]) <= a[i] ? out = true : out = out
    end

    return out

end

## MeshParameters are
## [yi,yf,lengthy,zi,zf,lengthz, gamma0,"x"]
function MakeMesh(SP::SystemParameters, MPs)

    gamma0     = MPs[7]
    alphaIndx  = range(MPs[1],stop = MPs[2], length = MPs[3])
    betaIndx   = range(MPs[4],stop = MPs[5], length = MPs[6])


    meshx =[]
    meshy =[]
    meshz =[]
    #This For defines the mesh within the supercell                                                                          \
                                                                                                                              
    for alphai in alphaIndx, betai in betaIndx
        
        if MPs[8] == "x"
            meshx = push!(meshx,gamma0)
            meshy = push!(meshy,alphai)
            meshz = push!(meshz,betai)
        elseif MPs[8] == "y"
            meshx = push!(meshx,alphai)
            meshy = push!(meshy,gamma0)
            meshz = push!(meshz,betai)
        elseif MPs[8] == "y"
            meshx = push!(meshx,alphai)
            meshy = push!(meshy,betai)
            meshz = push!(meshz,gamma0)
        end
            
            
    end

    return hcat(meshx,meshy,meshz)

end

function ECartesian(x,y,z,qlm,CenterSep)

    Etemp=[0.0 0.0 0.0]

    for i=1:2
    #Calculated in the Local Spherical Coordinate System
    # input x,y,z  for ESphrical in the global cartesian system
    Er     = ESpherical(x,y,z,i,qlm,CenterSep,m)[1]
    Etheta = ESpherical(x,y,z,i,qlm,CenterSep,m)[2]
    Ephi   = ESpherical(x,y,z,i,qlm,CenterSep,m)[3]

    # Transforms from Local Spherical Coordinate System to
    # global cartesian system.
    Etemp+=LSVFtoGCVF(i,CenterSep,Er,Etheta,Ephi,x,y,z)
    end

    return Etemp
end

#----------------------function------------------------------------
function Er(x,y,z,i,mm,nn,qlm,lmax,kvec,SP::SystemParameters)
    
    Rs            = SP.BaseVecs
    Er_at_rvec    = 0.0
    ts            = SP.LattVecs
    r0vec = (Rs[i,:]+mm*ts[1,:]+nn*ts[2,:])
  
    rvec  =  (Rs[i,:]+mm*ts[1,:]+nn*ts[2,:]) - [x,y,z] 
        
    R     = CartToSphere(rvec)[1]
    theta = CartToSphere(rvec)[2]
    phi   = CartToSphere(rvec)[3]

    QsPerParticle = Int(length(qlm)/SP.NP)
        
    for j = 1: QsPerParticle
        indx1 = (i-1)*QsPerParticle + j
        Ls = SP.IndxF(j,1,lmax,SP)
        l  = Ls[1]
        m  = Ls[3]
        qtemp =  qlm[indx1]*exp(im*dot(kvec,r0vec))   
        Er_at_rvec += InteractionMatrix.SphY(l,m,theta,phi)*qtemp*4*pi*(l+1)/( (2*l+1)*R^(l+2) )
            
    end
        
    return  Complex{Float64}(Er_at_rvec)
end

#----------------------function------------------------------------

function Eth(x,y,z,i,mm,nn,qlm,lmax,kvec,SP::SystemParameters)
    
    Rs            = SP.BaseVecs
    Eth_at_rvec    = 0.0
    ts             = SP.LattVecs

    r0vec = (Rs[i,:]+mm*ts[1,:]+nn*ts[2,:])
    rvec  =  (Rs[i,:]+mm*ts[1,:]+nn*ts[2,:]) - [x,y,z] 
        
    R     = CartToSphere(rvec)[1]
    theta = CartToSphere(rvec)[2]
    phi   = CartToSphere(rvec)[3]

       
    QsPerParticle = Int(length(qlm)/SP.NP)
        
    for j = 1: QsPerParticle
        indx1 = (i-1)*QsPerParticle + j
        Ls = SP.IndxF(j,1,lmax,SP)
        l  = Ls[1]
        m  = Ls[3]
        
        qtemp =  qlm[indx1]*exp(im*dot(kvec,r0vec))   
    
        m==1 && l==1 ? continue : m = m
            
        Eth_at_rvec += InteractionMatrix.DthSphY(l,m,theta,phi)*(-4.0)*qtemp/((2*l+1)*R^(l+2) )
            
    end
        
    
    return  Complex{Float64}(Eth_at_rvec)
end

#----------------------function------------------------------------

function Ephi(x,y,z,i,mm,nn,qlm,lmax,kvec,SP::SystemParameters)
    
    Rs            = SP.BaseVecs
    Ephi_at_rvec    = 0.0
    ts             = SP.LattVecs
    r0vec = (Rs[i,:]+mm*ts[1,:]+nn*ts[2,:])
   
    rvec  =  (Rs[i,:]+mm*ts[1,:]+nn*ts[2,:]) - [x,y,z] 
        
    R     = CartToSphere(rvec)[1]
    theta = CartToSphere(rvec)[2]
    phi   = CartToSphere(rvec)[3]

        
    QsPerParticle = Int(length(qlm)/SP.NP)
        
    for j = 1: QsPerParticle
        indx1 = (i-1)*QsPerParticle + j
        Ls = SP.IndxF(j,1,lmax,SP)
        l  = Ls[1]
        m  = Ls[3]
        
        qtemp =  qlm[indx1]*exp(im*dot(kvec,r0vec))   
    
        theta == 0 ? continue : m = m

        Ephi_at_rvec += -4pi*m*im*InteractionMatrix.SphY(l,m,theta,phi)*qtemp/( sin(theta)*(2*l+1)*R^(l+2) )
            
    end
    
    return  Complex{Float64}(Ephi_at_rvec)
end



#################################################
function ESpherical(x,y,z,i,mm,nn,qlm,lmax,k0vec,SP::SystemParameters)
    return [Er(x,y,z,i,mm,nn,qlm,lmax,k0vec,SP::SystemParameters) Eth(x,y,z,i,mm,nn,qlm,lmax,k0vec,SP::SystemParameters) Ephi(x,y,z,i,mm,nn,qlm,lmax,k0vec,SP::SystemParameters)]
end


function ECartesian2(x,y,z,qlm,lmax,polarization,SP::SystemParameters)

    Etemp=[0.0,0.0,0.0]
    NPs = SP.NP
    
    if InParticle([x,y,z],SP ::SystemParameters) == true 
        return [0.0 0.0 0.0] 
    end
    
    #Calculated in the Local Spherical Coordinate System
    # input x,y,z  in the global cartesian system
    for i = 1:NPs
        Es     = ESpherical(x,y,z,i,0,0,qlm,lmax,[0,0,0],SP::SystemParameters) 
        
        Er     = Es[1]
        Etheta = Es[2]
        Ephi   = Es[3]

        # Transforms from Local Spherical Coordinate System to
        # global cartesian system.
        Etemp += LSVFtoGCVF(i,Er,Etheta,Ephi,x,y,z,SP::SystemParameters)
    end
    
    E0   = SolveSystem.ExtPlaneWaveCart(polarization,SP::SystemParameters)

    return Etemp+E0
end

function ECartesian(x,y,z,qlm,NCells_EF,lmax,polarization,k0vec,SP::SystemParameters)

    Etemp=[0.0,0.0,0.0]
    NPs = SP.NP
    
    if InParticle_Ribbon([x,y,z],SP ::SystemParameters) == true 
        return 0.0
    end
    
    #Calculated in the Local Spherical Coordinate System
    # input x,y,z  in the global cartesian system
    for i = 1:NPs
         # Just for 1-D Ribbons
        for mm = -NCells_EF:NCells_EF, nn =  -NCells_EF:NCells_EF
            Es     = ElectricFieldMod.ESpherical(x,y,z,i,mm,nn,qlm,lmax,k0vec,SP::SystemParameters) 
        
            Er     = Es[1]
            Etheta = Es[2]
            Ephi   = Es[3]

            # Transforms from Local Spherical Coordinate System to
            # global cartesian system.
            Etemp += LSVFtoGCVF_periodic(i,mm,nn,Er,Etheta,Ephi,x,y,z,SP::SystemParameters)
        end
    end
    
    E0   = SolveSystem.ExtPlaneWaveCart(polarization,SP::SystemParameters)

#    return norm(Etemp+E0*exp(im*dot(k0vec,[x,y,z])))^2/ norm(E0*exp(im*dot(k0vec,[x,y,z])))^2
     return norm(Etemp)^2

end



function PlotEF(;omega_eV="200",polarization = "",lmax ="200",MP = "",SP::SystemParameters = "")
    
    ## MeshParameters are
    ## [yi,yf,lengthy,zi,zf,lengthz, x0]
    mesh = MakeMesh(SP::SystemParameters, MP)
    Qs   = SolveSystem.SolveQsAtLambda(lte(omega_eV),lmax,polarization,SP::SystemParameters)
    Ef(x,y,z) = norm(ECartesian(x,y,z,Qs,lmax,polarization,SP::SystemParameters) )^2
    
    EatMesh = map(Ef,mesh[:,1],mesh[:,2],mesh[:,3])    
    
    return hcat(mesh,EatMesh)
    
end
##--------------------------------------------------------------------
##    Electric Field Functions for Periodic Systems
##--------------------------------------------------------------------
function LSVFtoGCVF_periodic(i,mm,nn,Ar,Atheta,Aphi,x,y,z,SP::SystemParameters)
    
    Rs = SP.BaseVecs[i,:]
    t1 = SP.LattVecs[1,:]
    t2 = SP.LattVecs[2,:]
    
    x0 = Rs[1]+mm*t1[1]+nn*t2[1]
    y0 = Rs[2]+mm*t1[2]+nn*t2[2]
    z0 = Rs[3]+mm*t1[3]+nn*t2[3]
    z    = z-z0
    x    = x-x0
    y    = y-y0
   
    thet = ElectricFieldMod.CartToSphere([x y z])[2]
    phi  = ElectricFieldMod.CartToSphere([x y z])[3]

    Ax =Ar*sin(thet)*cos(phi) + Atheta*cos(thet)*cos(phi) - Aphi*sin(phi)
    Ay =Ar*sin(thet)*sin(phi) + Atheta*cos(thet)*sin(phi) + Aphi*cos(phi)
    Az =Ar*cos(thet)          - Atheta*sin(thet)

    Ax = Complex{Float64}(Ax)
    Ay = Complex{Float64}(Ay)
    Az = Complex{Float64}(Az)
   
    return [Ax, Ay, Az]
end



function InParticle_periodic(rvec,SP ::SystemParameters)
    t1 = SP.LattVecs[1,:]
    t2 = SP.LattVecs[2,:]
    t3 = t1+t2

    Rs = SP.BaseVecs
    a  = SP.Radius
    
    out = false 
    
    for mm = -3:3, nn = -3:3
        temp = mm*t1 + nn*t2
        for i = 1 : length(Rs[:,1])
            norm(rvec-Rs[i,:]-temp) <= a[i] ? out = true : out = out
            #norm(rvec-Rs[i,:]-t1) <= a[i] ? out = true : out = out
            #norm(rvec-Rs[i,:]-t2) <= a[i] ? out = true : out = out
            #norm(rvec-Rs[i,:]-t3) <= a[i] ? out = true : out = out

        end
    end
    
    return out
    
end

function HorizontalMesh(SP::SystemParameters, MeshParameters::Array )
    
    count   = 0
    t1 = SP.LattVecs[1,:]
    t2 = SP.LattVecs[2,:]
    t3 = t1+t2
   
    ll = Point(0,0)
    lr = Point(t1[1],t1[2])
    ur = Point(t3[1],t3[2])
    ul = Point(t2[1],t2[2])
    poly = Polygon(ll, lr, ur, ul)
    
    Xspan = range(-t1[1],stop = t1[1], length = Int(MeshParameters[1]) )
    Yspan = range(0,     stop = t3[2], length = Int(MeshParameters[1]) )

    Meshx = []
    Meshy = []
    for xi in Xspan, yi in Yspan

        InCell = inpolygon(poly, Point(xi,yi))
        InCell == true ? Meshx = push!(Meshx,xi) : xi = xi
        InCell == true ? Meshy = push!(Meshy,yi) : xi = xi
    
    end
    
    return hcat(Meshx,Meshy,fill(MeshParameters[3],length(Meshx)) )

end

function MakeMesh_SuperCell(SP::SystemParameters, MPs)

    z0  = MPs[2]
    ts  = SP.LattVecs

    Indx1   = range(0.01,stop = 1, length = MPs[1])
    Indx2   = range(0.01,stop = 1, length = MPs[1])


    meshx =[]
    meshy =[]
    meshz =[]
    #This For defines the mesh within the supercell                                                                          \
                                                                                                                              
    for mm in Indx1, nn in Indx2
        ri    = mm*ts[1,:]+nn*ts[2,:]
        
        meshx = push!(meshx,ri[1])
        meshy = push!(meshy,ri[2])
        meshz = push!(meshz,z0)
    end

    return hcat(meshx,meshy,meshz)

end

function MakeMesh_Ribbon(SP::SystemParameters, MPs)

    ts  = SP.LattVecs[1,:]
    Rs  = SP.BaseVecs
    yi  = minimum(Rs[:,2])
    yf  = maximum(Rs[:,2])    

    Indx_x   = range(0,  stop = norm(ts), length = MPs[1])
    Indx_y   = range(yi, stop = yf, length = MPs[2])
    z0  = MPs[3]


    meshx =[]
    meshy =[]
    meshz =[]
    #This For defines the mesh within the supercell                                                                          \
                                                                                                                              
    for mm in Indx1, nn in Indx2
        ri    = mm*ts[1,:]+nn*ts[2,:]
        
        meshx = push!(meshx,ri[1])
        meshy = push!(meshy,ri[2])
        meshz = push!(meshz,z0)
    end

    return hcat(meshx,meshy,meshz)

end

function ECartesian_periodic(x,y,z,qlm,NCells_EF,lmax,polarization,k0vec,SP::SystemParameters)

    Etemp=[0.0,0.0,0.0]
    NPs = SP.NP
    
    if InParticle_periodic([x,y,z],SP ::SystemParameters) == true 
        return 0.0
    end
    
    #Calculated in the Local Spherical Coordinate System
    # input x,y,z  in the global cartesian system
    for i = 1:NPs
        for mm = -NCells_EF:NCells_EF, nn =  -NCells_EF:NCells_EF
            Es     = ElectricFieldMod.ESpherical(x,y,z,i,mm,nn,qlm,lmax,k0vec,SP::SystemParameters) 
        
            Er     = Es[1]
            Etheta = Es[2]
            Ephi   = Es[3]

            # Transforms from Local Spherical Coordinate System to
            # global cartesian system.
            Etemp += LSVFtoGCVF_periodic(i,mm,nn,Er,Etheta,Ephi,x,y,z,SP::SystemParameters)
        end
    end
    
    E0   = SolveSystem.ExtPlaneWaveCart(polarization,SP::SystemParameters)

    return norm( Etemp+E0*exp(im*dot(k0vec,[x,y,z])) )
    #return norm(Etemp+E0)

end

function PlotEF_periodic(;omega_eV="",k0vec = "",polarization = "",NCells_EF= "",RepeatCell = "",lmax ="",MP = "",SP::SystemParameters = "")
    
    t1 = SP.LattVecs[1,:]
    t2 = SP.LattVecs[2,:]
    lambda = lte(omega_eV)
    

    ## MeshParameters are
    ## [Nx,Nywww,z0]
    mesh  = HorizontalMesh(SP, MP) #MakeMesh(SP, MPs)
    #mesh =  MakeMesh_SuperCell(SP::SystemParameters, MP)
    
    
    Qs   = SolveSystem.SolveQsAtLambda(lambda,lmax,polarization,NCells_EF,k0vec,SP::SystemParameters)

    Ef(x,y,z) = ECartesian_periodic(x,y,z,Qs,NCells_EF,lmax,polarization,k0vec,SP::SystemParameters)
    
    data    = map(Ef,mesh[:,1],mesh[:,2],mesh[:,3])    
    dataOut = data
    
    if RepeatCell == false
    
        return hcat(mesh,dataOut)
    
    else
        meshx = mesh[:,1]
        meshy = mesh[:,2]
        meshz = mesh[:,3]

        for mm = 0:1,nn = 0:1
            
            nn == 0 && mm == 0 ? continue : nn = nn
            
            rtemp = nn*t1+mm*t2
            xtemp = mesh[:,1].+rtemp[1]
            ytemp = mesh[:,2].+rtemp[2]
            
            dataOut  = vcat(dataOut,data)
            meshx    = vcat(meshx,xtemp)
            meshy    = vcat(meshy,ytemp)

        end
    
        return hcat(meshx,meshy,dataOut)
    
    end
end


function PlotEF_2Dperiodic(;omega_eV="",k0vec = "",polarization = "",NCells_EF= "",RepeatCell = "",lmax ="",MP = "",Nint = "",SP::SystemParameters = "")
    
    t1 = SP.LattVecs[1,:]
    t2 = SP.LattVecs[2,:]
    lambda = lte(omega_eV)

    Nint  = NCells_EF
    mesh  = MakeMesh_EF2D( MP,SP) 
    Qs    = SolveSystem.SolveQs_Periodic(lambda,lmax ,polarization ,k0vec ,Nint,SP::SystemParameters)

    Ef(x,y,z) = ElectricFieldMod.ECartesian_periodic(x,y,z,Qs,NCells_EF,lmax,polarization,k0vec,SP::SystemParameters)
    
    data    = map(Ef,mesh[:,1],mesh[:,2],mesh[:,3])    
    dataOut = data
    
    if RepeatCell == false
    
        return hcat(reshape(dataOut,Int(MP[3]),Int(MP[6]) ))
    
    else
        meshx = mesh[:,1]
        meshy = mesh[:,2]
        meshz = mesh[:,3]

        for mm = 0:1,nn = 0:1
            
            nn == 0 && mm == 0 ? continue : nn = nn
            
            rtemp = nn*t1+mm*t2
            xtemp = mesh[:,1].+rtemp[1]
            ytemp = mesh[:,2].+rtemp[2]
            
            dataOut  = vcat(dataOut,data)
            meshx    = vcat(meshx,xtemp)
            meshy    = vcat(meshy,ytemp)

        end
    
        return hcat(meshx,meshy,dataOut)
    
    end
end


function MakeMesh_EF2D( MPs,SP::SystemParameters)

    ts  = SP.LattVecs[1,:]+SP.LattVecs[2,:]
 

    Indx_x   = range(MPs[1],  stop = MPs[2] , length = Int(MPs[3]) )
    Indx_y   = range(MPs[4],  stop = MPs[5] , length = Int(MPs[6]) )
    
    z0  = MPs[7]
    
    meshx = []
    meshy = []
    meshz = []
    
    for i = 1:Int(MPs[3]), j = 1:Int(MPs[6])
        meshx = push!(meshx,Indx_x[i])
        meshy = push!(meshy,Indx_y[j])
        meshz = push!(meshz,z0)
    end

                                                                       
    if MPs[8] == "H"                                                                                                                        
        return hcat(meshx,meshy,meshz)
    else
       return  hcat(meshz,meshx,meshy)
    end
end

end