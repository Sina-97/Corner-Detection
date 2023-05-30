using Pkg   
# Pkg.add("ImageEdgeDetection")
# Pkg.add("Noise")
# Pkg.add("Images")
# Pkg.add("Plots")
# Pkg.add("ImageFiltering")
# Pkg.add("Random")
# Pkg.add("Statistics")
# Pkg.add("ImageEdgeDetection")
# Pkg.add("LinearAlgebra")
using Noise 
using Images
using Plots
using ImageFiltering
using Random
using Statistics
using ImageEdgeDetection
using ImageEdgeDetection: Percentile
using LinearAlgebra

plot_image(I; kws...) = plot(I; aspect_ratio=:equal, size=size(I), framestyle=:none, kws...)
img=load("blocks_bw.png")


function find_corners1(img,N,t,k)  
 Gx, Gy = imgradients(img, KernelFactors.sobel) 
 xi=[]
 yj=[]
func1=x->sum(x.^2)
func2=x->sum(x)
    window=(N,N)
Gx_Q=mapwindow(func1,Gx,window)
Gy_Q=mapwindow(func1,Gy,window)
Gxy_Q=mapwindow(func2,Gx.*Gy,window)
R=zeros(size(img,1),size(img,2))
R=((Gx_Q.*Gy_Q)-(Gxy_Q.^2))-k.*(Gx_Q+Gy_Q).^2
R_thresh=R.*(R.>t)
maxR=mapwindow(maximum,R_thresh,window)

for i=1:size(img,1)
    for j=1:size(img,2)
    if (R[i,j]==maxR[i,j]) & (R[i,j]>0)
        push!(xi,i)
        push!(yj,j)
    end
    end
end
return xi,yj
end

function find_corners2(img,N,t,k)  
    Gx, Gy = imgradients(img, KernelFactors.sobel) 
       xi=[]
       yj=[]
   func1=x->sum(x.^2)
   func2=x->sum(x)
       window=(N,N)
   Gx_Q=mapwindow(func1,Gx,window)
   Gy_Q=mapwindow(func1,Gy,window)
   Gxy_Q=mapwindow(func2,Gx.*Gy,window)
   eig=zeros(size(img,1),size(img,2))
   for i=1:size(img,1)
    for j=1:size(img,2)
        v=eigvals([Gx_Q[i,j] Gxy_Q[i,j];Gxy_Q[i,j] Gy_Q[i,j]])
        eig[i,j]=v[1]
    end
   end
   eig_thresh=eig.*(eig.>t)
   maxlanda=mapwindow(maximum,eig_thresh,window)
   for i=1:size(img,1)
    for j=1:size(img,2)
        if (eig_thresh[i,j]==maxlanda[i,j]) & (eig_thresh[i,j]>0)
            push!(xi,i)
            push!(yj,j)
        end
    end
   end
   return xi,yj
   end



x,y=find_corners1(img,5,0.02,0.04)
p1=plot(Gray.(img))
p1=plot!(y,x, seriestype=:scatter, label="corners")
display(p1)

x,y=find_corners2(img,5,0.065,0.04)
p2=plot(Gray.(img))
p2=plot!(y,x, seriestype=:scatter, label="corners")
display(p2)

# The result of these two methods are both good and they both were able to detect the corners,however, the threshold for these two 
# methods are different.
