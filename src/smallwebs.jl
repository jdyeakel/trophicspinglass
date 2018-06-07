function smallwebs(S)
    
    #Trophic chain
    Atc = zeros(Bool,S,S);
    Atc[diagind(Atc,-1)]=true;
    niche = collect(1:S);
    
    return Atc,niche
end
