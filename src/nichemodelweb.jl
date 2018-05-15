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
    keep = find(!iszero,sum(adjmatrix,2));
    
    while length(keep) < ladj
        ladj = size(adjmatrix)[1];
        keep = find(!iszero,sum(adjmatrix,2));
        adjmatrix = adjmatrix[keep,keep];
    end
    
    
    return adjmatrix
    
end

# R"""
# image($(adjmatrix),col=grey(c(0,1)))
# """
