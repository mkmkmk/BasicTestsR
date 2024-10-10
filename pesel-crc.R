setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())

# ad PESEL
w = c(1, 3, 7, 9, 1, 3, 7, 9, 1, 3, 1)
all = c("80051608853", "81013000106", "11111111116", "06660666002")
for (p in all)
{
    p = as.integer(unlist(strsplit(p, "")))
    print((10 - (sum(p[1:10]*w[1:10]) %% 10)) %% 10)
    stopifnot((10 - (sum(p*w) %% 10)) %% 10 == 0)
}
