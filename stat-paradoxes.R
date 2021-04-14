# https://scitechdaily.com/paradoxes-of-probability-statistical-strangeness/


################################################
# https://en.wikipedia.org/wiki/Simpson%27s_paradox#Vector_interpretation

# szkodliwe leczenie było zastosowane często dla najłatwierszych przypadków
# a rzadko dla najtrudniejszych przypadków
# i jednocześnie brak leczenia był zastosowany często dla trudnych przypadków 
# i rzadko dla łatwych przypadków
# dlatego wydało się że szkodliwe leczenie jest lepsze niż brak leczenia


a_num = c(10, 20, 30, 60)
a_rec = c(.1, .2, .3, .4)

a_q = a_num * a_rec
a_p = a_num

sum(a_num * a_rec) / 120
sum(a_num * a_rec)
a_num * a_rec


b_num = c(60, 30, 20, 10)
b_rec = c(.15, .3, .45, .6)
sum(b_num * b_rec) / 120
sum(b_num * b_rec)
b_num * b_rec

b_q = b_num * b_rec
b_p = b_num


# mimo że suma jest do 120 wykres do 60
plot(NA, xlim=c(0, 60), ylim=c(0,60))
for(i in 1:4)
    lines( c(0,a_p[i]), c(0,a_q[i]), col = 'red4')
lines(c(0, 120), c(0, sum(a_num * a_rec)), col = 'red4', lwd = 2, lty=3)


for(i in 1:4)
    lines( c(0,b_p[i]), c(0,b_q[i]), col = 'blue4')
lines(c(0, 120), c(0, sum(b_num * b_rec)), col = 'blue4', lwd = 2, lty=3)


###########################
# ad ratefallacy
300*.04
12*.08
(300-12)*.25

((300-12)*.25+12*.08)/300


###########################
# Multiple comparisons fallacy
# (Birthday Paradox)
1/365

sapply(2:50, function(n) n*(n-1)/2)

40*39/2 / 365

(1 - 1/365)^(40*39/2)
1 / (1 - 1/365)^(40*39/2)

(1-1/365)^(24*23/2)
(1-1/365)^(23*22/2)
(1-1/365)^(22*21/2)

sapply(2:50, function(n) n*(n-1)/2)

# p-stwo że w żadnej z tych wielu par nikt nie jednocześnie urodzin
pp = sapply(2:66, function(n) (1-1/365)^(n*(n-1)/2))
plot(2:66, pp, type = 'l')


