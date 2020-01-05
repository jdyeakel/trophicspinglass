loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/loadfuncs.jl");


#Build 

S=100;
C=0.01;
@time A,n = nichemodelweb(S,C);
#R"""
#image($(A),col=grey(c(0,1)))
#"""


S=50; C=0.1;
A,niche = smallwebs(S);
A,niche = nichemodelweb(S,C);





S=size(A)[1];
tl = trophic(A); tlsp = sortperm(tl,rev=true); tlsort = tl[tlsp];
Asort = A[tlsp,tlsp];
tmax=100;
#Predation coupling
kout=2.;
#Consumption coupling
kin=2.;
#Global influence of primary producers
gprim = 1.;
#Noise
sigma=0.01;
@time s = cascade(Asort,kout,kin,gprim,sigma,tmax);
# R"""
# image($(s),col=c('black','white'))
# """

#Change in state over time
#Low, middle high trophic levels
tlout = 8;
#tlcut = [round.(Int64,collect(1:(S/(tlout-1)):S));S];
#Delta =zeros(Int64,tlout-1,tmax-1);
#for i=1:tlout-1
#    tlrange1 = tlcut[i]; tlrange2 = tlcut[i+1];
#    for t=2:tmax
#        Delta[i,t-1] = sum((s[t-1,tlrange1:tlrange2] .- s[t,tlrange1:tlrange2]).^2);
#    end
#end


namespace = string("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/spin.pdf");
R"""
library(RColorBrewer)
library(igraph)
pal = brewer.pal(9,"Blues")
#pal = rev(brewer.pal($tlout-1,'Spectral'))
pdf($namespace,width=8,height=10)
par(mfrow=c(2,1))
image(x=$(collect(1:tmax)),y=$(collect(1:S)),$(s),col=c(pal[2],pal[7]),xlab='Time',ylab='Species ID')
#plot($(Delta[1,:]),type='l',ylim=c(min($Delta),max($Delta)),col=pal[1],lwd=2)
g = graph_from_adjacency_matrix($Asort)
plot(g,vertex.size=2,edge.arrow.size=0.2,vertex.label=NA,vertex.color='lightblue')
dev.off()
"""

#for i=tlout-1:-1:2
#    R"lines($(Delta[i,:]),col=pal[$i],lwd=2)"
#end
# R"dev.off()"
#Pattern of traveling peaks shows the cascade direction/strength






#t=50;
#plotweb(Asort,s[t,:],tlsort)






#webs
#tend = 100;
#trange = collect(tend-4+1:tend);
#R"par(mfrow=c(2,2))"
#for t=tend-4+1:tend
#    plotweb(Asort,s[t,:],tlsort)
#end
