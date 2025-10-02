if("rstudioapi" %in% rownames(installed.packages()) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


# Funkcja przeliczająca wysokość na ciśnienie według atmosfery standardowej
wysokosc_na_cisnienie <- function(H) {
    p0 <- 1013.25  # ciśnienie na poziomie morza w hPa
    p <- p0 * (1 - H/44300)^5.256
    return(p)
}

# Tworzenie danych
wysokosci <- seq(0, 20000, by = 100)  # wysokości od 0 do 20 km co 100m
cisnienia <- wysokosc_na_cisnienie(wysokosci)

# Tworzenie ramki danych
dane <- data.frame(
    wysokosc_m = wysokosci,
    wysokosc_km = wysokosci/1000,
    cisnienie_hPa = cisnienia
)

# Wykres
library(ggplot2)

p1 <- ggplot(dane, aes(x = wysokosc_km, y = cisnienie_hPa)) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = c(1013.25, 850, 700, 500, 300, 200),
               linetype = "dashed", alpha = 0.6, color = "red") +
    geom_vline(xintercept = c(1.5, 3, 5.5, 9, 12),
               linetype = "dashed", alpha = 0.6, color = "red") +
    labs(
        title = "Zależność ciśnienia od wysokości - Atmosfera Standardowa",
        subtitle = expression(paste("p = p"[0], " × (1 - H/44300)"^"5.256")),
        x = "Wysokość [km]",
        y = "Ciśnienie [hPa]",
        caption = "p₀ = 1013.25 hPa"
    ) +
    scale_x_continuous(breaks = seq(0, 20, 2)) +
    scale_y_continuous(breaks = c(0, 200, 300, 500, 700, 850, 1013.25)) +
    annotate("text", x = 16, y = 850, label = "850 hPa (~1.5 km)", size = 3) +
    annotate("text", x = 16, y = 700, label = "700 hPa (~3 km)", size = 3) +
    annotate("text", x = 16, y = 500, label = "500 hPa (~5.5 km)", size = 3) +
    annotate("text", x = 16, y = 300, label = "300 hPa (~9 km)", size = 3) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
    )

print(p1)

# Funkcja odwrotna - ciśnienie na wysokość
cisnienie_na_wysokosc <- function(p) {
    p0 <- 1013.25
    H <- 44300 * (1 - (p/p0)^(1/5.256))
    return(H)
}

# Przykładowe obliczenia dla poziomów meteorologicznych
poziomy_met <- c(1013.25, 850, 700, 500, 300, 200, 100)
wysokosci_met <- cisnienie_na_wysokosc(poziomy_met)

cat("Poziomy meteorologiczne:\n")
for(i in 1:length(poziomy_met)) {
    cat(sprintf("%7.2f hPa = %6.0f m = %5.2f km\n",
                poziomy_met[i], wysokosci_met[i], wysokosci_met[i]/1000))
}

# Wykres z zaznaczonymi poziomami meteorologicznymi
p2 <- ggplot(dane, aes(x = cisnienie_hPa, y = wysokosc_km)) +
    geom_line(color = "darkgreen", size = 1.2) +
    geom_point(data = data.frame(cisnienie = poziomy_met, wysokosc = wysokosci_met/1000),
               aes(x = cisnienie, y = wysokosc),
               color = "red", size = 3) +
    geom_text(data = data.frame(cisnienie = poziomy_met[2:6], wysokosc = wysokosci_met[2:6]/1000),
              aes(x = cisnienie, y = wysokosc,
                  label = paste0(cisnienie, " hPa")),
              hjust = -0.1, vjust = 0.5, size = 3) +
    labs(
        title = "Standardowe poziomy meteorologiczne",
        x = "Ciśnienie [hPa]",
        y = "Wysokość [km]"
    ) +
    scale_x_reverse() +  # odwrócona oś X (ciśnienie maleje z wysokością)
    theme_minimal()

print(p2)