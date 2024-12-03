if (F)
{
    rm(list = ls())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    print(getwd())
}

# t.test 2 groups

n = 500
group1 = rnorm(n / 2, 162.77, 20)
group2 = rnorm(n / 2, 168.11, 33)

dt = rbind(data.frame(y=group1, x = 1), data.frame(y=group2, x = 2))
dt = dt[order(dt$y),]

pop = par(mfrow=c(3,1))

plot(group1, col = 'red', type = 'p')
lines(group2, col = 'blue', type = 'p')

plot(1:n, dt$y, col=c('red','blue')[dt$x])

h1 = hist(group1, plot=F)
h2 = hist(group2, plot=F)
rgx = range(c(h1$mids, h2$mids)); rgy = c(-1, 1 + max(c(h1$counts, h2$counts)))
plot(h1$mids, h1$count, type='o', col='red', xlim = rgx, ylim = rgy)
lines(h2$mids, h2$count, type='o', col='blue')

par(pop)

print(t.test(y ~ x, data = dt))

shapiro.test(group1)
shapiro.test(group2)

var.test(group1, group2)

shapiro.test(rnorm(100, mean = 5, sd = 3))
shapiro.test(runif(100, min = 2, max = 4))
