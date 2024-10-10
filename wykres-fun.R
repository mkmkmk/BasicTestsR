setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())


x = seq(-100,100,0.2)
y = abs(x/2-1)
y1 = (x/2)

rx = c(-5,5)
ry = c(-4,4)

nx = diff(rx)
ny = diff(ry)

plot(x, y, type = 'p', col = 'blue',
     xlim = rx,
     ylim = ry,
     xaxp=c(rx, diff(rx)),
     yaxp=c(ry, diff(ry)),
)
abline(v = seq(rx[1], rx[2]), col = "gray", lty = "dotted")
abline(h = seq(ry[1], ry[2]), col = "gray", lty = "dotted")
lines(x, y1, type = 'p', col = 'red', xlim = rx, ylim=ry, xaxp=c(rx, diff(rx)), yaxp=c(ry, diff(ry)))


legend(x = -4, y = 4,
       legend = c("y = |x/2-1|", "y = x/2"),  # Nazwy funkcji
       col = c("blue", "red"),           # Kolory odpowiadające funkcjom
       pch = 1)






