#' ---
#' title: piasek jako magazyn energii
#' ---
#+ echo=F, include=F, warnings=F
# author: M. Krej
#-------------------

#' ### Magazyn energii 138 kWh, akumulator LiFePO4 PL
#' 
#' https://allegro.pl/oferta/magazyn-energii-138-kwh-akumulator-lifepo4-pl-11601879148
#' 
#' 379000 zł
#' 
#' 138 kWh
#' 
#' 
#' energia [MJ] :
138e3*3600/1e6 # MJ

#+ echo=F, include=F, warnings=F
#160Ah * 24V
160*24*3600/1e6 # MJ
496.8/13.8
138e3/160/24

#' ### piasek jako magazyn energii
#' ciepło właściwe : 
#' 
#' https://pl.wikipedia.org/wiki/Ciep%C5%82o_w%C5%82a%C5%9Bciwe
#' 
#' ciepło właściwe piasku [J/kg/K] :
c = 800 # J/kg/K
#' różnica temperatur [K] :
dT = 400 # K
#' energia [J] :
dQ = 500e6 # J

#' potrzebna masa [kg] :
dQ/c/dT # kg

#' ### sprawdzenie (w odwrotną stronę) : 
#' 
#' masa [kg] :
m = 1500 # kg

#' energia [MJ] :
c * m * dT / 1e6 # MJ

#' energia [kWh] :
c * m * dT /3600/1000 # kWh


#' czyli ten magazyn za 400k zł
#' magazynuje tyle energii co podgrzać 1.5 tony piasku o 400 stopni 
#' (nie licząc tego że to nie jest energia elektryczna)
#' 
#' ### +fiński magazyn energii: 
#' 
#' https://spidersweb.pl/2022/07/magazynowanie-energii-piasek.html
#' 



