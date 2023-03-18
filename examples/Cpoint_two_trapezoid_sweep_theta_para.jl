using SharedArrays
using Distributed
using DelimitedFiles
#addprocs(60)
@everywhere using RigorousCoupledWaveAnalysis

clear()
Air = ConstantPerm(1)
n2p02 = ConstantPerm(2.02^2)
a0 = 1000
N = 6
#θ = 1E-5 #elevation angle, zero will yield a singularity inversion error
α = 0 #azimuth angle
dx = 0.016
dy = 0
nx, ny, dnx, dny = ngrid(N, N)
alpha = 0.6
L0 = 247 / 450
L1 = sqrt(2 * L0^2 / (2 - alpha))
L2 = L1 * (1 - alpha)
H = L1
f = zeros(1000, 1000)
xi = range(-0.5, 0.5, 1000)
yi = range(-0.5, 0.5, 1000)
aC1 = 2 * H / (L2 - L1)
bC1 = H / 2 - H * L2 / (L2 - L1)
aC2 = -aC1
bC2 = bC1
for i = 1:length(xi)
    for j = 1:length(yi)
        x = xi[i]
        y = yi[j]
        #if -H/2 < y < H/2 && L2 * y + L1 * H > L1 * y + 2 * H * x && L2 * y + L1 * H > L1 * y - 2 * H * x
        if -H / 2 <= y <= H / 2 && y <= aC1 * x + bC1 && y <= aC2 * x + bC2
            f[i, j] = 1
        end
    end
end
#figure();
#imagesc(f);
F = real2recip(dnx, dny, f)
R1 = Custom(F)
R2 = Shift(R1, dx, dy)

h = 0.39 * a0
Layer1 = PatternedLayer(h / 2, [Air, n2p02], [R1])
Layer2 = PatternedLayer(h / 2, [Air, n2p02], [R2])
mdl = RCWAModel([Layer1, Layer2], Air, Air)

#f0 = 0.765:0.00005:0.784
f0 = 0.6:0.0005:0.9
wls = a0 ./ f0

#wls = 1050:1:1080
θi = -15+1E-5:0.05:15+1E-5
#θi = 1e-5
Rrf = SharedArray{Float64}(length(wls), length(θi))
#Rrf = zeros(length(wls), length(θi)) #Forward rcp reflectivity
Trf = SharedArray{Float64}(length(wls), length(θi))
Tlf = SharedArray{Float64}(length(wls), length(θi))
Rlf = SharedArray{Float64}(length(wls), length(θi))
inds = CartesianIndices(size(Rrf))

@time @sync @distributed for ind in inds
    #@time @inbounds Threads.@threads for ind in inds
    #@time for ind in inds
    i, j = ind.I
    λ = wls[i] #get wavelength from array
    θ = θi[j]
    grd = rcwagrid(N, N, a0, a0, θ, α, λ, Air) #build a reciprocal space grid
    ste, stm = rcwasource(grd, 1) #define source
    Rlf[i, j], Tlf[i, j] = RigorousCoupledWaveAnalysis.etm_reftra(sqrt(0.5) * (stm + 1im * ste), mdl, grd, λ) #lcp propagation
    Rrf[i, j], Trf[i, j] = RigorousCoupledWaveAnalysis.etm_reftra(sqrt(0.5) * (1im * stm + ste), mdl, grd, λ) #rcp propagation
end


#print(Tlf)
#print(Trf)

figure()
subplot(1, 2, 1)
imagesc(Tlf, clims = [0, 1], cmap = "gray")
#colorbar()
subplot(1, 2, 2)
imagesc(Trf, cmap = "gray", clims = [0, 1])
#colorbar()
writedlm("c_at_G_Trf.txt", Trf, '\t')
writedlm("c_at_G_Tlf.txt", Tlf, '\t')