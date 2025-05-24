

oblicz_sume_kontrolna <- function(numer) {
  numer <- gsub("/", "", toupper(numer))
  
  if (nchar(numer) != 12) {
    stop("Numer KW (bez cyfry kontrolnej) musi mieć dokładnie 12 znaków.")
  }
  
  wartosci_liter <- c(
    X = 10, A = 11, B = 12, C = 13, D = 14, E = 15,
    F = 16, G = 17, H = 18, I = 19, J = 20, K = 21,
    L = 22, M = 23, N = 24, O = 25, P = 26, R = 27,
    S = 28, T = 29, U = 30, W = 31, Y = 32, Z = 33
  )
  
  znaki <- strsplit(numer, "")[[1]]
  
  wartosci <- sapply(znaki, function(z) {
    if (grepl("[0-9]", z)) {
      as.numeric(z)
    } else if (z %in% names(wartosci_liter)) {
      wartosci_liter[[z]]
    } else {
      stop(paste("Nieprawidłowy znak w numerze KW:", z))
    }
  })
  
  wagi <- rep(c(1, 3, 7), 4)
  suma <- sum(wartosci * wagi)
  suma <- suma %% 10
  return(suma)
}


oblicz_sume_kontrolna("WA1M/00000003") # /5
