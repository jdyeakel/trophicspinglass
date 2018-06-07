function trophic(A)
    R"""
    library(MASS)  
    library(NetIndices)
    rtl<-TrophInd($(A'))
    """
    @rget rtl;
    tl = Array(rtl[1]) - 1;
    return tl
end

    