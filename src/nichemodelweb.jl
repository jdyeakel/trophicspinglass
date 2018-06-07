function nichemodelweb(S,C) 
    
    
    #list of random niche values
    n = sort(rand(S));
    
    #Define beta distribution
    a = 1;
    b = (1/(2*C)) - 1;
    bdist = Beta(a,b);
    
    #range
    r = n.*rand(bdist,S);
    
    #center of range
    c = rand.(Uniform.(r/2,n));
    
    #Find species that fall within the range of each other
    prey = Array{Array{Int64}}(S);
    [prey[i] = find(x-> (x > c[i] - (r[i]/2) && x < c[i] + (r[i]/2)), n) for i=1:S];
    

    adjmatrix = zeros(Bool,S,S);
    [adjmatrix[i,prey[i]] = true for i=1:S];
    
    ladj = size(adjmatrix)[1];
    keep = find(!iszero,vec(sum(adjmatrix,2))+vec(sum(adjmatrix,1)));
    niche = n;
    
    
    while length(keep) < ladj
        ladj = size(adjmatrix)[1];
        keep = find(!iszero,vec(sum(adjmatrix,2))+vec(sum(adjmatrix,1)));
        adjmatrix = adjmatrix[keep,keep];
        niche = niche[keep];
    end
    
    
    
    return adjmatrix', niche
    
end

# R"""
# image($(adjmatrix),col=grey(c(0,1)))
# """
