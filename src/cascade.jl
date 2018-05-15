function cascade(S,C,kout,kin,sigma,tmax)
    
    A = nichemodelweb(S,C);
    n = size(A)[1];
    
    s = Array{Int64}(tmax,n);
    #Define initial state
    s0 = rand([-1,1],n);
    s[1,:] = s0;

    G = DiGraph(A');
    
    for t=2:tmax
        past_s = copy(s[t-1,:]);
        future_s = zeros(Int64,n);
        #Determine the futures state of node i
        for i=1:n
            
            #Consumers of species i
            #Consumer influence is negative
            nn_out = outneighbors(G,i);
            
            #Resources of species i
            #Resource influence is positive
            nn_in = inneighbors(G,i);
            
            #Spin function
            newstate = sign((-kout*sum(past_s[nn_out]) + kin*sum(past_s[nn_in])) + sigma*rand(epsilondist));
            
            future_s[i] = newstate;
        end
        s[t,:] = future_s;
    end
    
    return(
    A',
    s
    )
    
end