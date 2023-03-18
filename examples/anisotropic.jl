using RigorousCoupledWaveAnalysis
using LinearAlgebra
no=1.52 #ordinary index
ne=1.81 #extraordinary index
#LC=ConstantPermA([no^2,1,1,ne^2,no^2]) #sparse permittivity tensor (exx,exy,eyx,eyy,ezz)
LC = ConstantPermA([ε_MoO3[1,1],ε_MoO3[1,2],ε_MoO3[2,1],ε_MoO3[2,2],ε_MoO3[3,3]])
#for other LC orientation, you would have to apply a rotation matrix to the first four elements
air=ConstantPerm(1.0)
#wls=500:5:1500 #wavelength array
omega_cm = 10:100:20000
#omega_cm = 10
wls = 1 ./ (omega_cm/1e-2) * 1e9

N=2
pitch=500
l1=AnisotropicLayer(1000,LC) #only 1um LC layer in air
mdl=RCWAModel([l1],air,air)
Rx=zeros(length(wls))
Ry=zeros(length(wls))
Tx=zeros(length(wls))
Ty=zeros(length(wls))
rx=zeros(length(wls))*1im
ry=zeros(length(wls))*1im
tx=zeros(length(wls))*1im
ty=zeros(length(wls))*1im
for i=1:length(wls) 
    grd=rcwagrid(N,N,500,500,1e-5,0,wls[i],air) #incoming wave is at 45 degrees polarization
    ste,stm=rcwasource(grd)
    kzin=grd.k0[3]#z component of incoming wave momentum
    ems=RigorousCoupledWaveAnalysis.eigenmodes(grd,wls[i],mdl.layers)#layer eigenmodes
    ref=RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,mdl.εsup,wls[i])#superstrate eigenmodes
    tra=RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,mdl.εsub,wls[i])#superstrate eigenmodes
    ro,to,r,t=etm_propagate(ref,tra,ems,stm,grd,false)#perform etm computation
    #need to create masks to select the 0-th order x and y
    maskx=zeros(length(ro))
    maskx[N*(2N+1)+N+1]=1
    masky=zeros(length(ro))
    masky[(2N+1)^2+N*(2N+1)+N+1]=1
    #compute powers
    Tx[i]=-RigorousCoupledWaveAnalysis.a2p(to.*maskx,to*0,tra.V,I,kzin) #power in 0-th order x transmission
    Ty[i]=-RigorousCoupledWaveAnalysis.a2p(to.*masky,to*0,tra.V,I,kzin) #power in 0-th order y transmission
    Rx[i]=RigorousCoupledWaveAnalysis.a2p(0*ro,ro.*maskx,ref.V,I,kzin) #power in 0-th order x reflection
    Ry[i]=RigorousCoupledWaveAnalysis.a2p(0*ro,ro.*masky,ref.V,I,kzin) #power in 0-th order y reflection                    
    #for insight into these check src/ETM/ETM.jl
    #compute amplitudes (only valid in lossless media)
    ex,ey=RigorousCoupledWaveAnalysis.a2e2d(to,I) #identity matrix of correct size I             
    #println(length(ex))
    tx[i]=ex[N*(2N+1)+N+1]
    ty[i]=ey[N*(2N+1)+N+1]
    ex,ey=RigorousCoupledWaveAnalysis.a2e2d(ro,I) #identity matrix of correct size
    rx[i]=ex[N*(2N+1)+N+1]
    ry[i]=ey[N*(2N+1)+N+1]
    #print(kzin)
end
#now we can also compute polarization angle (the real part of this is the angle, the imaginary part is related to ellipticity:
phi_R=atand.(ry./rx).-45 #rotation in reflected wave
phi_T=atand.(ty./tx).-45 #rotation in transmitted wave
figure()
plot(wls,Rx)
plot(wls,Ry)