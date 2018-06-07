function cascade(A,kout,kin,gprim,sigma,tmax)
    
    # A,niche = nichemodelweb(S,C);
    n = size(A)[1];
    
    s = Array{Int64}(tmax,n);
    #Define initial state
    s0 = rand([-1,1],n);
    s[1,:] = s0;

    G = DiGraph(A);
    
    epsilondist = Normal(0,1);
    
    g = 0;
    
    for t=2:tmax
        
        past_s = copy(s[t-1,:]);
        future_s = zeros(Int64,n);
        
        #Determine the futures state of node i
        for i=1:n
            
            #PREDATION
            #Predator influence is negative
            nn_out = outneighbors(G,i);
            
            #Resource consumption
            #Resource influence is positive
            nn_in = inneighbors(G,i);
            
            if length(nn_in) == 0
                g = gprim;
            end
            
            #Spin function
            newstate = sign(-kout*sum(past_s[nn_out]) + kin*sum(past_s[nn_in]) + sigma*rand(epsilondist) + g);
            
            # newstate = sign(-kin*sum(past_s[nn_in]) + sigma*rand(epsilondist));
            
            future_s[i] = newstate;
        end
        
        # nn_out = outneighbors.(G,collect(1:n));
        # nn_in = inneighbors.(G,collect(1:n));
        # past_out = Array{Array}(n);
        # past_in = Array{Array}(n);
        # eprand = rand(epsilondist,n);
        # for i=1:n
        #     past_out[i] = past_s[nn_out[i]];
        #     past_in[i] = past_s[nn_in[i]];
        # end
        # #Spin function
        # newstate = sign.(-kout.*sum.(past_out) + kin.*sum.(past_in) + sigma.*eprand);
        # future_s = newstate;
        
        s[t,:] = future_s;
    end
    
    return(
    s
    )
    
end
