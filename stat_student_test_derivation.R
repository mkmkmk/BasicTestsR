if (F)
{
    rm(list = ls())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    print(getwd())
}


# 1. For a single sample mean, we know that:
#     X̄ ~ N(μ, σ²/n)
#     → mean has normal distribution with the same mean and n-times smaller variance
#
# 2. For difference of means (under H₀: μ₁ = μ₂):
#     (X̄₁ - X̄₂) ~ N(0, σ₁²/n₁ + σ₂²/n₂)
#     → variance of differences equals: σ₁²/n₁ + σ₂²/n₂
#
# 3. Standardizing (dividing by standard error):
#     Z = (X̄₁ - X̄₂) / √(σ₁²/n₁ + σ₂²/n₂) ~ N(0,1)
#
# 4. We don't know σ₁ and σ₂, so we use variance estimators:
#     s² = Σ(Xᵢ - X̄)² / (n-1)
#     → such variable has chi-square distribution,
#       specifically (n-1)s²/σ² (here σ²=1) has chi-square distribution with (n-1) degrees of freedom
#
# 5. Substituting estimators:
#     t = (X̄₁ - X̄₂) / √(s₁²/n₁ + s₂²/n₂)
#     → such random variable has Student's t-distribution
#       because a normally distributed variable divided by square root of a
#       chi-square distributed variable (divided by degrees of freedom)
#       has t-distribution


# -------------------
# ad 1. X̄ ~ N(μ, σ²/n)
{
    n = 30
    pop_mean = 10
    pop_sd = 2
    sample_means = replicate(100000, mean(rnorm(n, pop_mean, pop_sd)))
    hist(sample_means, freq=FALSE, breaks = 100)
    curve(dnorm(x, pop_mean, pop_sd/sqrt(n)), add=TRUE, col="blue", lwd=2)
}


# -------------------
# ad 2. (X̄₁ - X̄₂) ~ N(0, σ₁²/n₁ + σ₂²/n₂)
{
    n = 30
    diff_means = replicate(100000, mean(rnorm(n)) - mean(rnorm(n)))
    hist(diff_means, freq=FALSE, breaks = 100)
    curve(dnorm(x, 0, sqrt(1/n + 1/n)), add=TRUE, col="blue", lwd=2)
}


# -------------------
# ad 4. (n-1)s²/σ² has chi-square distribution with (n-1) degrees of freedom
{
    n = 8
    chi_sq_stats = replicate( 100000, {
        x = rnorm(n)  # ~ N(0,1)
        (n-1)*var(x)/1  # (n-1)s²/σ², σ²=1
    })
    hist(chi_sq_stats, freq=FALSE, breaks=100)
    curve(dchisq(x, df=n-1), add=TRUE, col="blue", lwd=2)
}


# -------------------
# ad 5. (X̄₁ - X̄₂) / √(s₁²/n₁ + s₂²/n₂ has Student's t-distribution
{
    n = 8
    t_stats = replicate(100000, {
        x = rnorm(n)
        mean(x)/(sd(x)/sqrt(n))
    })
    hist(t_stats, freq=FALSE, breaks = 100, xlim = c(-10,10))
    curve(dt(x, df=n-1), add=TRUE, col="blue", lwd=2)
}


# -------------------
# ad 5b. normally distributed variable Z divided by square root of a
# chi-square distributed variable V (divided by degrees of freedom)
# has t-distribution t = Z/√(V/(n-1))
{
    n = 1000000
    df = 5
    Z = rnorm(n)
    V = rchisq(n, df)
    t = Z/sqrt(V/df)
    hist(t, freq=FALSE, breaks=250, xlim=c(-10,10))
    curve(dt(x, df=df), add=TRUE, col="blue", lwd=2)
}


# -------------------
# implementation of the two-group t-test
{
    my_ttest = function(x1, x2) {

        n1 = length(x1)
        n2 = length(x2)

        mean1 = mean(x1)
        mean2 = mean(x2)

        var1 = var(x1)
        var2 = var(x2)

        t_stat = (mean1 - mean2) / sqrt(var1/n1 + var2/n2)

        # degrees of freedom (Welch)
        df = (var1/n1 + var2/n2)^2 /
            ((var1^2/(n1^2*(n1-1))) + (var2^2/(n2^2*(n2-1))))

        p_val = 2*pt(-abs(t_stat), df=df)

        return(list(t_stat=t_stat, df=df, p_val=p_val))
    }

    x1 = rnorm(30, 10, 2)
    x2 = rnorm(30, 12, 2)

    print(t.test(x1, x2))
    print(my_ttest(x1, x2))
}


# -------------------
# with increasing degrees of freedom, the Student distribution becomes a normal distribution
{
    curve(dnorm(x), -4, 4, col="black", lwd=2)
    curve(dt(x, df=1), add=TRUE, col="red", lty=2)    # Cauchy'ego
    curve(dt(x, df=2), add=TRUE, col="blue", lty=2)
    curve(dt(x, df=5), add=TRUE, col="green", lty=2)
    curve(dt(x, df=30), add=TRUE, col="purple", lty=2)
}


# -------------------
# visualization of the distribution tails
{
    df = 10
    t_stat = 2.5

    curve(dt(x, df), from=-4, to=4)
    abline(v=c(-t_stat, t_stat), col="red", lty=2)
    abline(h=0)

    x = seq(-4, -abs(t_stat), length=100)
    polygon(c(x, rev(x)), c(dt(x, df), rep(0, length(x))), col="lightblue")
    x = seq(abs(t_stat), 4, length=100)
    polygon(c(x, rev(x)), c(dt(x, df), rep(0, length(x))), col="lightblue")

    # p-value
    p_value = 2 * pt(-abs(t_stat), df)
    print(paste("p-value =", p_value))
}


# -------------------
# a sample from a normal distribution has a normal distribution
{
    big = rnorm(1000000)
    sub = sample(big, 10000)

    h1  = hist(big, plot = F)
    h2 = hist(sub, plot = F)

    rgx = range(c(h1$mids, h2$mids));
    plot(h1$mids, h1$density, type='o', col='red', xlim = rgx)
    lines(h2$mids, h2$density, type='o', col='blue')

    sd(big)
    sd(sub)
}
