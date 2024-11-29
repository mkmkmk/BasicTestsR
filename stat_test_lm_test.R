if (F)
{
    rm(list = ls())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    print(getwd())
}

x = seq(1, 5, 0.1)
y = 5e2 * x + 3e3 + runif(length(x))*5e3

model = lm(y ~ x)
plot(x, y, ylim = range(y) + 1.5 * c(-1,1) * diff(range(y)), col = 'blue')
abline(model, col = "red")  # dodanie linii trendu

cor.test(x, y)

print(summary(model))
print(anova(model))

print(anova(model)$"F value"[1])
print(summary(model)$fstatistic[1])