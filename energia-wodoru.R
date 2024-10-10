#' ---
#' title: ad energia wodoru
#' ---
#' 

#' udział masy wodoru w wodzie
# stosunek mas wodoru i tlenu
# H * 2 + O
# masa atomowa wodoru:1 masa tlenu: 16
# 1 * 2 + 16 --> 1 + 8
#' udział tlenu:
8 / (8 + 1)
#' udział wodoru:
1 / (8 + 1)

#' 1 litr wody (1000g) zawiera tyle gramów wodoru
1000 * 1 / (8 + 1)

#' "One litre of hydrogen weighs 0.09 gram", ale na Wikipedii:
#' 
#' https://pl.wikipedia.org/wiki/Wod%C3%B3r
#' 
#' 0.082 kg/m3 == g/litr
#'
#' z jednego litra wody (1000g) tyle litrów wodoru
1000 * 1 / (8 + 1) / 0.09

#' z jednego grama wody tyle litrów wodoru:
1 / (8 + 1) / 0.09


#' Wh energii zużyto:
7.5 * 0.06 * 1

#' J energii zużyto:
7.5 * 0.06 * 3600

#' gram wody przerobionej na gaz
0.36 

#' "It is possible to produce 1.23 litres of hydrogen from one gram of water"
#' 
#' uzyskano litrów wodoru: 
0.36 * 1.23

#' gęstość energii wodoru MJ/litr: 0.01079
#' 
#' wytworzony wodór ma tyle J energii:
0.36 * 1.23 * 0.01079 * 1e6

#' sprawność %
(0.36 * 1.23 * 0.01079 * 1e6) / (7.5 * 0.06 * 3600) * 100
(0.36 * 1 / (8 + 1) / 0.082 * 0.01079 * 1e6) / (7.5 * 0.06 * 3600) * 100


#+ echo=F, include=F, warnings=F
#brudnopis:
10.1/0.01079
1000/0.09



