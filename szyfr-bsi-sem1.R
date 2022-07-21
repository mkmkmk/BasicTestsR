if (FALSE)
  install.packages("combinat")

if (FALSE)
  rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())

library(combinat)

inp = gsub(" ", "", "╣ ⌠ ┼ ﴿  ╡ ╬ ╪ ├ ╠  ╬ Ł ⌠ ╓ ⌠ Ł  ╡ ┬ ┌ ╠ ═ ╖ ─ ┼  ┬ ╖ ═ ╓ ┌ ╡")

a = c("A", "S")
b = c("C", "F", "I", "T", "Y")
c = c("E", "H", "K", "L", "M", "O", "R", "Z")

all = expand.grid(permn(a), permn(b), permn(c))
nrow(all)

a = c("⌠", "╡")
b = c("╣", "﴿", "╪", "─", "├")
c = c("┼", "╬", "╠", "╓", "┬", "┌", "═", "╖")

#fileConn = file("output.txt", "w")

# ii=1
# for(ii in 1:100)
for(ii in 1:nrow(all))
{
  row = all[ii,]
  
  A = row[[1]][[1]]
  B = row[[2]][[1]]
  C = row[[3]][[1]]
  
  ares = inp
  # i = 1
  for(i in 1:length(A))
    ares = gsub(a[i], A[i], ares)
  for(i in 1:length(B))
    ares = gsub(b[i], B[i], ares)
  for(i in 1:length(C))
    ares = gsub(c[i], C[i], ares)
  
  #cat(sprintf("%s\n",ares), file = fileConn, append = TRUE)
  if ((length(grep("taki", ares, ignore.case = T)) > 0) && (length(grep("szyfr", ares, ignore.case = T)) > 0))
      cat(sprintf("%10d: %s\n", ii, ares))
      
  if ((ii %% 100e3) == 0)
      cat(sprintf("%10d\n", ii))
}

#close(fileConn)

