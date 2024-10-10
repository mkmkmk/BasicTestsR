nm = c("1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p")
sm = c(2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6)
sum = 0
names = c()
for(i in 1:length(nm))
{
    names = c(names, nm[i])
    sum = sum + sm[i]
    cat(sprintf("%3d: %s\n", sum, paste(names, collapse = ", ")))
}

