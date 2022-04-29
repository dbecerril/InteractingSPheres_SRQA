module InitializeSystem

using DielectricFunctions
using Plots
using LinearAlgebra
using SystemCoordinates

export SystemParameters,DefineSquareRibbon,DefineRibbonSystem
export DefineFiniteSystem,Define2dPeriodicSystem,DefineRotatedBilayerSystem
export DefineZigZagRibbon,DefineZigZagRibbon_dy,DefineBZigZagRibbon
export DefineArmChairRibbon,DefineBArmChairRibbon,DefineSC
export DimerMap,ViewBilayerSystem,ViewSystem,DefineAABilayer,RR
export NParticleMap,MatIndxToL,MatIndxToM,MatIndxToI,RotateVectorList
export AA_base,AB_base,TriangLVecs,clip_zrange

mutable struct SystemParameters
    
    Radius     ::Array{Float64}                     
    Emed       ::Float64
    Specie     ::String
    BaseVecs   ::Array{Float64}
    LattVecs   ::Array{Float64}
    NP         ::Int
    IndxF      ::Function
    LayerSep   ::Float64
    SystemName ::String
    NNSep      ::Float64
    RotIndx    ::Array{Float64}
    
end

RR(Theta) = [cos(Theta) -sin(Theta) 0 ; sin(Theta) cos(Theta) 0; 0 0 1]

#----------------------------------------------------------------------------
#---- Set of Functions to initialize Different Types of Systems
#----------------------------------------------------------------------------
function DefineFiniteSystem(;SystemName::String = "", radii = "",BaseVecs = "",Specie ="",Emed ="",IndxF::Function = "" )
    NP = length(radii)
    
    if length(radii) != length(BaseVecs[:,1])
        error("Number of Particles not equal to number of Base Vectors \n")
    else
        return SystemParameters(radii,Emed,Specie,BaseVecs,[0 0 0;0 0 0],NP,IndxF,0.0,SystemName,1e5,[0,0])
    end
end

function Define2dPeriodicSystem(;SystemName::String= "", radii = "",BaseVecs = "",LattVecs = "",Specie ="",Emed ="",IndxF::Function = "",NNSep = "")
    NP = length(radii)
    
    if length(radii) != length(BaseVecs[:,1])
        error("Number of Particles not equal to number of Base Vectors \n")
    else
        return SystemParameters(radii,Emed,Specie,BaseVecs,LattVecs,NP,IndxF,0.0,SystemName,NNSep,[0,0])
    end
end

function DefineSC(;SystemName::String= "", radius = "",Specie ="",Emed ="",IndxF::Function = "",NNSep = "", SCP = "")
    a      = radius
    NSC    = Int(SCP[1])
    aa     = fill(a,2*NSC^2)
    
    theta1 = SCP[2]
    vv     = sqrt(3)*(NNSep)*[0.5 sqrt(3)/2 0.0;-0.5 sqrt(3)/2 0.0  ]
    t1     = RR(theta1)*vv[1,:]
    t2     = RR(theta1)*vv[2,:]

    vv1    = vcat(transpose(t1),transpose(t2))
    d      = RR(theta1)*[0.0,NNSep,0.0]

    bb = []

    for m = 0:NSC-1, n = 0:NSC-1
        bA = n*t1+m*t2
        bB = n*t1 + m*t2 + d
        bb = append!(bb,bA)
        bb = append!(bb,bB)
    end

    bb = transpose(reshape(bb,3,:))

    MonolayerSC = Define2dPeriodicSystem(;SystemName= SystemName, radii =aa ,BaseVecs = bb,LattVecs = NSC*vv1,Specie =Specie,Emed =Emed,IndxF = NParticleMap,NNSep = NNSep )

    return MonolayerSC
end


function DefineAABilayer(;SystemName::String= "", radii = "",BaseVecs = "",LattVecs = "",LayerSep = "",Specie ="",Emed ="",IndxF::Function = "", NNSep = "")
    
    NP = length(radii)
    if length(radii) != length(BaseVecs[:,1])
        error("Number of Particles not equal to number of Base Vectors \n")
    else
        return SystemParameters(radii,Emed,Specie,BaseVecs,LattVecs,NP,IndxF,LayerSep,SystemName,NNSep,[0,0])
    end
end

function DefineRotatedBilayerSystem(;SystemName = "",Radius = "",dgap = "" ,LayerSep = "", Specie ="",Emed ="",IndxF::Function = "",RotIndx="",NNSep = "")
        
    geometry = RotatedBilyerCoordinates(Radius,dgap,RotIndx,LayerSep);
    
    Top = geometry[1]
    Bot = geometry[2]
    L1  = geometry[3]
    L2  = geometry[4]
    NP  = length(Top[:,1])+length(Bot[:,1])

    return SystemParameters(fill(Radius,NP),Emed,Specie,vcat(Top,Bot),vcat(transpose(L1),transpose(L2) ),NP,IndxF,LayerSep,SystemName,NNSep,RotIndx)
    
end

function DefineZigZagRibbon(;SystemName = "", Radius = "",CSep = "" ,Specie ="",Emed ="",NWidth = "" )

    BVecs = ZigZagRibbonCoordinates(CSep,Radius,NWidth);
    aa    = fill(Radius,length(BVecs[:,1]))
    LVecs = [norm(sqrt(3)*CSep) 0 0;0 0 0]
    SPout =  Define2dPeriodicSystem(; radii = aa,BaseVecs = BVecs,LattVecs = LVecs, Specie =Specie,Emed =Emed,IndxF = NParticleMap,SystemName = SystemName,NNSep = CSep )

    return SPout
    
end

function DefineZigZagRibbon_dy(;SystemName = "", Radius = "",CSep = "" ,Specie ="",Emed ="",NWidth = "",dy = "" )

    BVecs = ZigZagRibbonCoordinates_dy(CSep,Radius,NWidth,dy);
    aa    = fill(Radius,length(BVecs[:,1]))
    LVecs = [norm(sqrt(3)*CSep) 0 0;0 0 0]
    SPout =  Define2dPeriodicSystem(; radii = aa,BaseVecs = BVecs,LattVecs = LVecs, Specie =Specie,Emed =Emed,IndxF = NParticleMap,SystemName = SystemName,NNSep = CSep+dy )

    return SPout
    
end

function DefineBZigZagRibbon(;SystemName = "", Radius = "",CSep = "" ,Specie ="",Emed ="",NWidth = "" )
    
    BVecs = BZigZagRibbonCoordinates(CSep,Radius,NWidth);
    aa    = fill(Radius,length(BVecs[:,1]))
    LVecs = [norm(sqrt(3)*CSep) 0 0;0 0 0]
    SPout =  Define2dPeriodicSystem(; radii = aa,BaseVecs = BVecs,LattVecs = LVecs, Specie =Specie,Emed =Emed,IndxF = NParticleMap,SystemName= SystemName,NNSep = CSep )

    return SPout
    
end

function DefineArmChairRibbon(;SystemName = "", Radius = "",CSep = "", Specie ="",Emed ="",NWidth = "" )
    
    BVecs = DefineArmChairCoordinates(CSep,Radius,NWidth)
    L1    = sqrt(3)*(CSep)*([0.5 sqrt(3)/2 0.0] + [-0.5 sqrt(3)/2 0.0  ])

    aa    = fill(Radius,length(BVecs[:,1]))
    LVecs = [norm(L1) 0 0;0 0 0]
    SPout =  Define2dPeriodicSystem(; radii = aa,BaseVecs = BVecs,LattVecs = LVecs, Specie =Specie,Emed =Emed,IndxF = NParticleMap, SystemName= SystemName,NNSep = CSep)

    return SPout
    
end

function DefineBArmChairRibbon(;SystemName ="", Radius = "",CSep = "", Specie ="",Emed ="",NWidth = "" )
    
    BVecs = DefineBArmChairCoordinates(CSep,Radius,NWidth)
    L1    = sqrt(3)*(CSep)*([0.5 sqrt(3)/2 0.0] + [-0.5 sqrt(3)/2 0.0  ])

    aa    = fill(Radius,length(BVecs[:,1]))
    LVecs = [norm(L1) 0 0;0 0 0]
    SPout =  Define2dPeriodicSystem(; radii = aa,BaseVecs = BVecs,LattVecs = LVecs, Specie =Specie,Emed =Emed,IndxF = NParticleMap,SystemName= SystemName,NNSep = CSep)

    return SPout
    
end

function DefineSquareRibbon(;SystemName::String ="", Radius = "",CSep = "", Specie ="",Emed ="",NWidth = "" )
    
    BVecs = [0 0 0]
    for i = 1 : NWidth - 1
        
        BVecs = vcat(BVecs,[0 i*CSep 0])
        
    end
    
    aa    = fill(Radius,NWidth)
    LVecs = [CSep 0 0;0 0 0]
    SPout =  Define2dPeriodicSystem(; radii = aa,BaseVecs = BVecs,LattVecs = LVecs, Specie =Specie,Emed =Emed,IndxF = NParticleMap,SystemName = SystemName,NNSep = CSep)

    return SPout
    
end



function DefineRibbonSystem(; SystemName = "", RibbonType = "", Radius ="" ,CSep = "" ,Specie = "",Emed = "", NWidth = "")
    if RibbonType == "ZZ"
        return DefineZigZagRibbon(;SystemName = SystemName, Radius = Radius ,CSep = CSep,Specie = Specie,Emed = Emed, NWidth = NWidth )
    elseif RibbonType == "BZZ"
        return DefineBZigZagRibbon(;SystemName = SystemName, Radius = Radius ,CSep = CSep,Specie = Specie,Emed = Emed, NWidth = NWidth )
    elseif RibbonType == "AC"
        return DefineArmChairRibbon(;SystemName = SystemName, Radius =Radius ,CSep = CSep,Specie = Specie,Emed = Emed, NWidth = NWidth )
    elseif RibbonType == "BAC"
        return DefineBArmChairRibbon(;SystemName = SystemName, Radius = Radius ,CSep = CSep,Specie = Specie,Emed = Emed, NWidth = NWidth )        
    elseif RibbonType == "SQ"
        return DefineSquareRibbon(;SystemName = SystemName, Radius = Radius ,CSep =CSep, Specie =Specie,Emed =Emed,NWidth = NWidth )
    else
        error("RibbonType must be: \n ZZ:zigzag \n BZZ:Bearded Zigzag \n AC: Armchair \n BAC:Bearded Armchair \n SQ: Square \n")
    end
end
#################################
### Index Maps for A dimer    ###
#################################
function DimerMap(Indx1::Int,Indx2::Int,lmax::Int)
    
    if Indx1 == 0 && Indx2 == 0  
        return 2*lmax
    else
        return DimerMatMap(Indx1,Indx2,lmax)    
    end
    
end

function DimerMatMap(Indx1,Indx2,lmax)
    #lfunctions
    l1 = ifelse(Indx1 <= lmax, Indx1, Indx1-lmax)
    l2 = ifelse(Indx2 <= lmax, Indx2, Indx2-lmax)
    
    #msfunctions
    m1 = 0
    m2 = 0
    
    #ifunctions
    i = ifelse(Indx1 <= lmax, 1, 2)
    j = ifelse(Indx2 <= lmax, 1, 2)
    
    return [l1,l2,m1,m2,i,j]
end

#####################################
### Index Maps for N Particles    ###
#####################################

function NParticleMap(Indx1::Int,Indx2::Int,lmax::Int,SP::SystemParameters)
    
    NPs = SP.NP
    
    if Indx1 == 0 && Indx2 == 0  
        return NPs*(lmax*(lmax+1)+lmax)
    else
        l1 = MatIndxToL(Indx1,lmax,SP)
        l2 = MatIndxToL(Indx2,lmax,SP)
        m1 = MatIndxToM(Indx1,lmax,SP)
        m2 = MatIndxToM(Indx2,lmax,SP)
        i  = MatIndxToI(Indx1,lmax,SP)
        j  = MatIndxToI(Indx2,lmax,SP)
        
        return [l1,l2,m1,m2,i,j]
    end
    
end

function MatIndxToL(MIndx,lmax,SP)
    ##################################################
    ### Locates in which SubMatrix the index Belongs to
    ##################################################
    SubMatDim      = lmax*(lmax+1)+lmax
    NPs = SP.NP
    
    for i = 1:NPs
        (i-1)*SubMatDim < MIndx <= i* SubMatDim ? MIndx = MIndx - (i-1)*SubMatDim : SubMatDim = SubMatDim 
        (i-1)*SubMatDim < MIndx <= i* SubMatDim ? break : SubMatDim = SubMatDim 
    end

    ##################################################
    ### Locates  which "L" it belongs to 
    ##################################################

    list1 = Int32[l*(l+1)+l for l = 1:lmax]
    out = 0
    
    lmax == 1 ? lmaxt = lmax+1 : lmaxt = lmax
    
    for i = 2:lmaxt
        if MIndx <= 3 
            out = 1
            break
        elseif list1[i-1] < MIndx <= list1[i]
            out = i 
            break
        end
        
    end
        
    return out
end


##########################################################
### Aux. Function later used to define H elements
### Input : "global" matrix index
### Output: "local" Angular Moment Index m
### SubMatDim variable refers to the dimension of a submatrix
### total Matrix Dim will be 2*MatDim
#########################################################
function MatIndxToM(MIndx,lmax,SP)
    
    ##################################################
    ### Locates in which SubMatrix the index Belongs to
    ##################################################
    SubMatDim      = lmax*(lmax+1)+lmax
    NPs = SP.NP
    
    for i = 1:NPs
        (i-1)*SubMatDim < MIndx <= i* SubMatDim ? MIndx = MIndx - (i-1)*SubMatDim : SubMatDim = SubMatDim 
        (i-1)*SubMatDim < MIndx <= i* SubMatDim ? break : SubMatDim = SubMatDim 
    end
    ##################################################
    ### Locates  which Ltemp "m" belongs in 
    ##################################################
    list1 = Int32[l*(l+1)+l for l = 1:lmax]
    out = 0
    lmax == 1 ? lmaxt = 2 : lmaxt = lmax
    
    for i = 2:lmaxt
        if MIndx <= 3 
            out =  1
            break
        elseif list1[i-1] < MIndx <= list1[i]
            out = i
            break
        end
        
    end

    l1 = out-1
    TempIndx = MIndx - ( l1*(l1+1)+l1 )
    #@show TempIndx out MIndx
    list2 = Int32[m for m = -out:out]
    
    return list2[TempIndx]
    
end


##########################################################
### Aux. Function later used to define H elements
### Input : "global" matrix index
### Output: "local" Angular Moment Index m
### SubMatDim variable refers to the dimension of a submatrix
### total Matrix Dim will be 2*MatDim
#########################################################
function MatIndxToI(MIndx,lmax,SP)
    
    ##################################################
    ### Locates in which SubMatrix the index Belongs to
    ##################################################
    SubMatDim      = lmax*(lmax+1)+lmax
    NPs = SP.NP
    out = 0
    
    for i = 1:NPs
        (i-1)*SubMatDim < MIndx <= i* SubMatDim ? out = i : SubMatDim = SubMatDim 
        (i-1)*SubMatDim < MIndx <= i* SubMatDim ? break : SubMatDim = SubMatDim 
    end
    
    return out
    
end


function CircleShape(x0,y0,r)
    theta = range(0,stop = 2pi,length = 30)
    
    return x0.+r*sin.(theta),y0.+r*cos.(theta)
    
end


function ViewSystem(SP::SystemParameters,NCells)    
    
    Rs  = SP.BaseVecs
    NPs = SP.NP
    as  = SP.Radius
    t1  = SP.LattVecs[1,:]
    t2  = SP.LattVecs[2,:]
    t3  = t1+t2
    Out = plot(CircleShape(Rs[1,1],Rs[1,2],as[1]),seriestype = :shape,aspect_ratio = 1,label = "",color = :blue,fmt = :png) 
    SuperCell = [0 0;t1[1] t1[2];t3[1] t3[2] ;t2[1] t2[2]; 0 0]

    # ------------- View 1-D Systems
    if norm(t2) == 0
        
        for nn = 1:NCells
        
            Rtemp = Rs[1,1:2]+nn*t1
            plot!(CircleShape(Rtemp[1],Rtemp[2],as[1]),seriestype = :shape,aspect_ratio = 1, c= :blue,label = "")
        end
    
        if NPs > 1
            for i = 2:NPs,nn = 0:NCells
        
                Rtemp = Rs[i,1:2]+nn*t1
                plot!(CircleShape(Rtemp[1],Rtemp[2],as[i]),seriestype = :shape,aspect_ratio = 1, c= :blue,label = "")

            end
        end
     
        return Out
        
    else

    for i = 1:NPs,nn = -NCells:NCells,mm = -NCells:NCells
        
        Rtemp = Rs[i,:]+nn*t1+mm*t2
        if Rs[i,3] == 0    
            plot!(CircleShape(Rtemp[1],Rtemp[2],as[i]),seriestype = :shape,aspect_ratio = 1, c= :blue,label = "")
        else
           plot!(CircleShape(Rtemp[1],Rtemp[2],as[i]),seriestype = :shape,aspect_ratio = 1, c= :red,alpha = 0.5,label = "")
        end
     end
         plot!(SuperCell[:,1],SuperCell[:,2],lab = "Super Cell",lw = 3, c = :green)

        return Out
        
    end
end



function RotateVectorList(VectorList,Theta)
    RR(Theta) = [cos(Theta) -sin(Theta) 0 ; sin(Theta) cos(Theta) 0; 0 0 1]
    
    NVecs = length(VectorList[:,1])
    Out   = zeros(NVecs,3) 
    
    for i = 1: NVecs
        Out[i,:] = transpose(RR(Theta)*VectorList[i,:])
    end
        
    return Out
end

function TriangLVecs(CSep)
    return sqrt(3)*(CSep)*[0.5 sqrt(3)/2 0.0;-0.5 sqrt(3)/2 0.0  ]
end

function AA_base(CSep,LSep)
    vv    = TriangLVecs(CSep)
    
    return [0 0 0;0 CSep 0;0 0 -LSep;0 CSep -LSep]
end

function AB_base(CSep,LSep)
    vv    = TriangLVecs(CSep)
    return [0 0 0;0 CSep 0;0 0 -LSep;0 -CSep+2*vv[1,2] -LSep]
end

function clip_zrange(z,zmin,zmax)
    
    if zmin == false && zmax == false
        return z
    else
        if z < zmin 
            return zmin
        elseif z>zmax
            return zmax
        else
            return z
        end
    end
    
end

end