

x = seq(-5, 5, length.out =  1000)
dist = exp(-x^2 / 2) / sqrt(2 * pi)
plot(x, dist, type = 'l')
runif

gaussMonte = function(n)
{
    sapply(1:n, 
            function(tmp)
            {
                # dla x na brzegach b. rzadko będą losowane y pod krzywą
                # dla x na środku prawie zawsze będą losowane y pod krzywą
                repeat
                {
                    x = runif(1, -10, 10)
                    y = runif(1, 0, 1 / sqrt(2 * pi))
                    if(y < exp(-x^2 / 2) / sqrt(2 * pi))
                        break
                }
                x        
            })   
}


gaussMonte(40)

if(F)
{
    hist(rnorm(10000), freq = F, breaks = 100)
    lines(x, dist, type = 'l', lwd=2)
    sd(rnorm(10000))
    shapiro.test(rnorm(10))
}

{
    to_ver = gaussMonte(10000)
    hist(to_ver, freq = F, breaks = 100)
    lines(x, dist, type = 'l', lwd=2)
    shapiro.test(to_ver[1:100])
    sd(to_ver)
    mean(to_ver)
}
