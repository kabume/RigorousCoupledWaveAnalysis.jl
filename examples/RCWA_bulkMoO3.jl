using RigorousCoupledWaveAnalysis
using PyPlot
using LinearAlgebra

#The parameters of MoO3
c = 299792458*1e2 #cm/s
hp = 6.62607e-034 #Js
ω = 925.9 #cm^-1
wl = 1 ./ ω*1e7
f = c .* ω

εx_inf = 5.78
ωx_LO1 = 534.3; ωx_TO1 = 506.7; Gammax1 = 49.1 
ωx_LO2 = 974.5; ωx_TO2 = 821.4; Gammax2 = 6.8 
ωx_LO3 = 999.2; ωx_TO3 = 998.7; Gammax3 = 0.35 
εx=εx_inf*(1+(ωx_LO1^2-ωx_TO1^2)./(ωx_TO1^2-ω.^2-im*ω*Gammax1)+ (ωx_LO2^2-ωx_TO2^2)./(ωx_TO2^2-ω.^2-im*ω*Gammax2)+(ωx_LO3^2-ωx_TO3^2)./(ωx_TO3^2-ω.^2-im*ω*Gammax3))

εy_inf = 6.07
ωy_LO = 850.1; ωy_TO = 544.6; Gammay = 9.5 
εy = εy_inf * (1+(ωy_LO^2-ωy_TO^2) ./ (ωy_TO^2-ω.^2-im*ω*Gammay))

εz_inf = 4.47
ωz_LO = 1006.9; ωz_TO=956.7; Gammaz=0.65
εz = εz_inf * (1+(ωz_LO^2-ωz_TO^2) ./ (ωz_TO^2-ω.^2-im*ω*Gammaz))
ε_MoO3 = [εx 0 0;0 εy 0;0 0 εz]
Δθ = 30
R = [cosd(Δθ) -sind(Δθ) 0; sind(Δθ) cosd(Δθ) 0; 0 0 1]
ε_MoO3 = transpose(R) * ε_MoO3 * R


#clear()
Air = ConstantPerm(1)
MoO3 = ConstantPermA([ε_MoO3[1,1],ε_MoO3[1,2],ε_MoO3[2,1],ε_MoO3[2,2],ε_MoO3[3,3]]) #(exx,exy,eyx,eyy,ezz)
a0 = 1 #nm
N = 1
θ = 1E-5 #elevation angle, zero will yield a singularity inversion error
α = 0 #azimuth angle

h = 300 * a0
Layer1 = AnisotropicLayer(h,MoO3)
mdl = RCWAModel([Layer1], Air, Air)

omega_cm = 10:100:20000
#omega_cm = 10
wls = 1 ./ (omega_cm/1e-2) * 1e9

r_xx = zeros(length(wls))*im
r_yy = zeros(length(wls))*im
r_xy = zeros(length(wls))*im
r_yx = zeros(length(wls))*im

#TM:Ex, TE, Ey
for i = eachindex(wls)
    λ = wls[i] #get wavelength from array
    grd = rcwagrid(N, N, a0, a0, θ, α, λ, Air) #build a reciprocal space grid
    ste, stm = rcwasource(grd) #define source
    ems = RigorousCoupledWaveAnalysis.eigenmodes(grd,λ,mdl.layers) #eigenmodes of all layers
    sup = RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,mdl.εsup,λ) #superstrate and substrate
    sub = RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,mdl.εsub,λ)
    kzin = grd.k0[3] #real(sqrt(get_permittivity(m.εsup,λ)))
    rtm,to,r,t=etm_propagate(sup,sub,ems,stm,grd,false) #propagate amplitudes
    rte,to,r,t=etm_propagate(sup,sub,ems,ste,grd,false) #propagate amplitudes

    exx,exy=RigorousCoupledWaveAnalysis.a2e2d(rtm,I)
    eyx,eyy=RigorousCoupledWaveAnalysis.a2e2d(rte,I)

    r_xx[i] = sum(exx)
    r_yy[i] = sum(eyy)
    r_xy[i] = sum(exy)
    r_yx[i] = sum(eyx)
#=
    maskx=zeros(length(rtm))
    maskx[N*(2N+1)+N+1]=1
    masky=zeros(length(rtm))
    masky[(2N+1)^2+N*(2N+1)+N+1]=1
    #compute powers
    r_xx[i]=RigorousCoupledWaveAnalysis.a2p(0*rtm,rtm.*maskx,sup.V,I,kzin) #power in 0-th order x reflection
    r_yx[i]=RigorousCoupledWaveAnalysis.a2p(0*rtm,rtm.*masky,sup.V,I,kzin) #power in 0-th order y 
=#
end

figure()
plot(omega_cm,abs.(r_xx).^2,"k-")
#plot(omega_cm,abs.(r_yy).^2,"r-")
#plot(omega_cm,abs.(r_xy).^2,"g--")
plot(omega_cm,abs.(r_yx).^2,"b--")
#plot(omega_cm,abs.(r_y).^2,'bo')

