#' ---
#' title: Wyprowadzenie wymiarów formatu arkusza A4 (w języku *R*)
#' author: M.K.
#' output: pdf_document
#' ---
#+ echo=F, include=F, warnings=F

#' ### Warunek proporcji
#'
#' Takie same proporcje boków po złożeniu kartki.
#'
#' Lewa strona: stosunek boku dłuższego $a$ do krótszego $b$ kartki złożonej na pół
#'
#' Prawa strona: tak samo, ale po **rozłożeniu** kartki, stosunek boku dłuższego $2b$ do krótszego $a$:
#'
#' $$\frac{a}{b} = \frac{2b}{a}$$
#'
#' ### Przekształcenia:
#'
#' $$a = \frac{2b^2}{a}$$
#'
#' $$a^2 = 2b^2$$
#'
#' ### Wynik:
#'
#' $$a = \sqrt{2} \cdot b$$
#'
#' ### Format A0
#'
#' Główne założenie serii "A" - powierzchnia równa $1 \text{ m}^2$:
#'
#' $$a \cdot b = 1$$
#'
#' $$\frac{a}{b} = \sqrt{2}$$
#'
#' $$a = \frac{1}{b}$$
#'
#' $$\frac{1}{b^2} = \sqrt{2}$$
#'
#' $$b^2 = \frac{1}{\sqrt{2}}$$
#'
#' $$b = \frac{1}{\sqrt[4]{2}}$$

b = 1 / sqrt(sqrt(2))
a = 1 / b

a
b

#' ### Funkcja realizująca składanie kartki na pół:
fold = function(x){ a = x[1]; b = x[2]; return(c(b, a/2)); }

#' ### Wymiary A4 [m]
fold(fold(fold(fold(c(a, b)))))

#' ### Wymiary A4 [mm]
as.integer(1000*fold(fold(fold(fold(c(a, b))))))
