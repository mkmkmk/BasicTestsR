

if(F)
{
    getSD = function(amp)
    {
        sig = amp * runif(1000, -1, 1)
        # plot(sig, type = 'l')
        sd(sig)
    }

    amp = 1:400
    resSD = sapply(amp, function(a) getSD(a))
    
    plot(amp, resSD)
    plot(amp,resSD/amp)
    mean(resSD / amp) 
    lines(amp, rep(0.58, length(amp)))
}


{
    amp = 33
    sig = amp * runif(10000, -1, 1)
    # plot(sig)
    
    amp *.58
    
    sd(sig)/amp
    
    sd(sig / (sd(sig) / .33))
    
    sd(sig / (amp *.58 / .33))
    # plot(sig / (amp *.58 / .33))
    
    quantize = function(samp)
    {
        if(samp > .5) 
            2 
        else if (samp < -.5) 
            -2 
        else if (samp > 0)
            1
        else 
            -1
    }
    
    hist(sapply(sig / (amp * .58 / .33), quantize))

}



