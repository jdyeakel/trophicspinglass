loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/loadfuncs.jl");


#Build 

S=1000;
C=0.1;
@time A,n = nichemodelweb(S,C);
R"""
image($(A),col=grey(c(0,1)))
"""


S=500;
C=0.02;
# A,niche = smallwebs(S);
A,niche = nichemodelweb(S,C);
S=size(A)[1];
tl = trophic(A); tlsp = sortperm(tl,rev=true); tlsort = tl[tlsp];
Asort = A[tlsp,tlsp];
tmax=500;
#Predation coupling
kout=0.0000;
#Consumption coupling
kin=1;
#Global influence of primary producers
gprim = 0.000;
#Noise
sigma=0.001;
@time s = cascade(Asort,kout,kin,gprim,sigma,tmax);
# R"""
# image($(s),col=c('black','white'))
# """

#Change in state over time
#Low, middle high trophic levels
tlout = 8;
tlcut = [convert.(Int64,round.(collect(1:(S/(tlout-1)):S),0));S];
Delta =zeros(Int64,tlout-1,tmax-1);
for i=1:tlout-1
    tlrange1 = tlcut[i]; tlrange2 = tlcut[i+1];
    for t=2:tmax
        Delta[i,t-1] = sum((s[t-1,tlrange1:tlrange2] .- s[t,tlrange1:tlrange2]).^2);
    end
end
R"""
library(RColorBrewer)
pal = rev(brewer.pal($tlout-1,'Spectral'))
par(mfrow=c(1,2))
image(x=$(collect(1:tmax)),y=$(collect(1:S)),$(s),col=c('black','white'))
plot($(Delta[1,:]),type='l',ylim=c(min($Delta),max($Delta)),col=pal[1],lwd=2)
"""
for i=tlout-1:-1:2
    R"lines($(Delta[i,:]),col=pal[$i],lwd=2)"
end
#Pattern of traveling peaks shows the cascade direction/strength






t=50;
plotweb(Asort,s[t,:],tlsort)






#webs
tend = 100;
trange = collect(tend-4+1:tend);
R"par(mfrow=c(2,2))"
for t=tend-4+1:tend
    plotweb(Asort,s[t,:],tlsort)
end
