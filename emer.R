if("rstudioapi" %in% rownames(installed.packages()) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


a = trunc(as.integer(as.Date('2025-01-30') - Sys.Date()) / 365.25 * (252 - 36))  
b = trunc(as.integer(as.Date('2025-01-30') - Sys.Date())) # 25.46

c = trunc(as.integer(as.Date('2028-02-15') - Sys.Date()) / 365.25 * (252 - 36) - 3 * 5)
d = trunc(as.integer(as.Date('2028-02-15') - Sys.Date()))

e = trunc(as.integer(as.Date('2024-08-16') - Sys.Date()) / 365.25 * (252 - 36))
f = trunc(as.integer(as.Date('2024-08-16') - Sys.Date()))

print(c(a, b, c, d, e, f))

