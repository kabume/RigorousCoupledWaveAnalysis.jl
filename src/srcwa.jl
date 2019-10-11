module srcwa
using LinearAlgebra
export srcwa_reftra,a2p,slicehalf,scatterSource,srcwa_matrices,Srcwa_matrices,srcwa_amplitudes,srcwa_abs
using ..models
using ..materials
using ..eigenmodes
using ..grid
using ..scatterMatrices

struct Srcwa_matrices
    Kzref::AbstractArray{Complex{Float64},2}
    Kztra::AbstractArray{Complex{Float64},2}
    Slayers::Array{ScatterMatrix,1}
end

function a2p(a,Kx,Ky,Kz,kz0)
    ex,ey,ez=a2e(a,Kx,Ky,Kz)
    return e2p(ex,ey,ez,Kz,kz0)
end
function a2e(a,Kx,Ky,Kz)
    ex,ey=slicehalf(a)
    ez=-Kz\(Kx*ex+Ky*ey)
    return ex,ey,ez
end
function a2e2d(a,W)
    e=W*a
    ex,ey=slicehalf(e)
    return ex,ey
end

function scatterSource(kinc,Nx,Ny)
    #the total number of scattering states
    width=(Nx*2+1)*(Ny*2+1)
    #vertical
    normal=[0,0,1]
    #te polarization E-field is perpendicular with z-axis and propagation direction (so, parallel with surface)
    kte=cross(normal,kinc)/norm(cross(normal,kinc))
    #tm polarization E-field is perpendicular with te and propagation direction (so, not necessarily parallel with surface)
    ktm=cross(kinc,kte)/norm(cross(kinc,kte))
    esource=zeros(width*2)*1im
    esource[convert(Int64,(width+1)/2)]=kte[1]
    esource[convert(Int64,(width+1)/2)+width]=kte[2]
    a0te=esource
    esource=zeros(width*2)*1im
    esource[convert(Int64,(width+1)/2)]=ktm[1]
    esource[convert(Int64,(width+1)/2)+width]=ktm[2]
    a0tm=esource#/sqrt(epsref)
    return a0te,a0tm
end
#just slices a vector e in half
function slicehalf(e)
    mylength=convert(Int64,size(e,1)/2)
    return e[1:mylength,:],e[mylength+1:end,:]
end
function e2p(ex,ey,ez,Kz,kz0)
    P=abs.(ex).^2+abs.(ey).^2+abs.(ez).^2
    P=sum(real.(Kz)*P/real(kz0))
    return P
end


function srcwa_matrices(model::Model,grd::Rcwagrid,λ)
    Slayers=Array{ScatterMatrix,1}(undef,length(model.layers)+2)
    for ct=2:length(Slayers)-1
        Slayers[ct]=scattermatrix_layer(eigenmode(grd.dnx,grd.dny,grd.Kx,grd.Ky,grd.k0,λ,model.layers[ct-1]),grd.V0)
    end
    ref=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsup,λ))
    tra=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsub,λ))
    Slayers[1]=scattermatrix_ref(ref,grd.V0)
    Slayers[end]=scattermatrix_tra(tra,grd.V0)
    return Srcwa_matrices(ref.Kz,tra.Kz,Slayers)
end
function srcwa_reftra(a0,grd::Rcwagrid,mtr::Srcwa_matrices)
    S=concatenate(mtr.Slayers)
    #a0te,a0tm=srcwa_source(grd.kin,Nx,Ny)
    aRte=S.S11*a0
    R=a2p(S.S11*a0,grd.Kx,grd.Ky,mtr.Kzref,grd.kin[3])
    T=a2p(S.S21*a0,grd.Kx,grd.Ky,mtr.Kztra,grd.kin[3])
    return R,T
end

function srcwa_amplitudes(a0,grd::Rcwagrid,mtr::Srcwa_matrices)
    a=zeros(length(a0),length(mtr.Slayers)-1)*1im
    b=zeros(length(a0),length(mtr.Slayers)-1)*1im
    for ct=1:size(a,2)
        Sbefore=concatenate(mtr.Slayers[1:ct])
        Safter=concatenate(mtr.Slayers[ct+1:end])
        a[:,ct]=(I-Sbefore.S22*Safter.S11)\(Sbefore.S21*a0)
        b[:,ct]=(I-Safter.S11*Sbefore.S22)\(Safter.S11*Sbefore.S21*a0)
    end
    return a,b
end

function srcwa_abs(a,b,grd::Rcwagrid)
    p=zeros(size(a,2))
    for ct=1:size(a,2)
        ex,ey=a2e2d(a[:,ct]+b[:,ct],I)
        hx,hy=a2e2d(-a[:,ct]+b[:,ct],grd.V0)
        p[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/grd.kin[3]
    end
    return p
end



end
