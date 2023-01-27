if("rstudioapi" %in% rownames(installed.packages()) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


sol = as.Date('2022-12-22')
day = as.Date('2023-01-12')
day = Sys.Date()

if (sol > day) sol + (sol - day) else sol - (day - sol)

for (dt in 0:51)
{
    day = sol - dt
    res = if (sol > day) sol + (sol - day) else sol - (day - sol)
    print(c(res, day))
}


