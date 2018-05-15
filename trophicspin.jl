using Distributions
using RCall
using LightGraphs

@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/nichemodelweb.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/cascade.jl")


#Build 

S=1000;
C=0.1;
@time A = nichemodelweb(S,C);
R"""
image($(A),col=grey(c(0,1)))
"""


S=1000;
C=0.1;
tmax=200;
kout=1;
kin=2;
sigma=0.1;
@time A,s = cascade(S,C,kout,kin,sigma,tmax);
R"""
image($(s),col=grey(c(0,1)))
"""


