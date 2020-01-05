loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/loadfuncs.jl");

gvec = collect(0:0.5:50);
lgvec = length(gvec);
koutvec = collect(0:0.5:15);
lkoutvec = length(koutvec);
reps = 200;
propcycle = SharedArray{Float64}(lgvec,lkoutvec,reps);

paramvec = Array{Int64}(undef,lkoutvec*lgvec*reps,3);
paramvec[:,1] = repeat(collect(1:lgvec),inner=lkoutvec*reps);
paramvec[:,2] = repeat(collect(1:lkoutvec),inner=reps,outer=lgvec);
paramvec[:,3] = repeat(collect(1:reps),outer=lkoutvec*lgvec);
tlvec = lgvec*lkoutvec*reps;

#Consumption coupling
kin=2.;

@time @sync @distributed for ii=1:tlvec
    
    i = paramvec[ii,1];
    j = paramvec[ii,2];
    r = paramvec[ii,3];
        
    gprim = gvec[i];
    kout = koutvec[j];

    S=100; C=0.01;
    # A,niche = smallwebs(S);
    A,niche = nichemodelweb(S,C);
    S=size(A)[1];
    tl = trophic(A); tlsp = sortperm(tl,rev=true); tlsort = tl[tlsp];
    Asort = A[tlsp,tlsp];
    tmax=500;
    #Predation coupling
    # kout=2.;
    #Consumption coupling
    # kin=2.;
    #Global influence of primary producers
    # gprim = 20.;
    #Noise
    sigma=0.01;
    s = cascade(Asort,kout,kin,gprim,sigma,tmax);
    
    #analysis
    
    #proportion of oscillating nodes
    burnin = 400;
    sdiff = vec(sum(diff(s[burnin:tmax,:],dims=1).!=0,dims=1));
    propcycle[i,j,r] = length(findall(x->x>0,sdiff))/size(s,2);
    
end
mpropcycle = mean(propcycle,dims=3)[:,:,1];

namespace = string("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/propcycle.pdf");
R"""
library(RColorBrewer)
library(fields)
pal=brewer.pal(9,"Spectral");
pdf($namespace,width=8,height=7)
image.plot(x=$gvec,y=$koutvec/$kin,z=$mpropcycle,col=pal,xlab='Primary productivity',ylab='Top-down influence')
dev.off()
"""



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
