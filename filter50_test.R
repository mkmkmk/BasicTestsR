rm(list = ls())


# Implementacja filtru w R
Filter50 <- function(fs) {
    Bw <- 9.95
    F50 <- 50

    theta <- 2 * pi * F50 / fs
    beta <- 2 * pi * Bw / fs
    alpha1 <- 2 * cos(theta) / (1 + tan(beta / 2))
    alpha2 <- (1 - tan(beta / 2)) / (1 + tan(beta / 2))

    num <- c((1 + alpha2) / 2, -alpha1, (1 + alpha2) / 2)
    den <- c(1, -alpha1, alpha2)

    list(num = num, den = den, fs = fs)
}

frequency_response <- function(filter, frequencies) {
    num <- filter$num
    den <- filter$den
    fs <- filter$fs

    omega <- 2 * pi * frequencies / fs

    # e^jω
    z <- exp(1i * omega)

    # Licznik: b0 + b1*z^(-1) + b2*z^(-2)
    numerator <- num[1] + num[2] * z^(-1) + num[3] * z^(-2)

    # Mianownik: 1 + a1*z^(-1) + a2*z^(-2)
    denominator <- den[1] + den[2] * z^(-1) + den[3] * z^(-2)

    H <- numerator / denominator

    magnitude <- abs(H)
    phase <- Arg(H) * 180 / pi

    data.frame(
        frequency = frequencies,
        magnitude = magnitude,
        magnitude_dB = 20 * log10(magnitude),
        phase = phase
    )
}


library(ggplot2)

# częstotliwość próbkowania 1 kHz
fs <- 1000
filter <- Filter50(fs)

frequencies <- seq(1, fs/2, length.out = 1000)

response <- frequency_response(filter, frequencies)
filter_tf <- list(b = filter$num, a = filter$den)
class(filter_tf) <- "Arma"

fqz = freqz(filter_tf, Fs = fs, n=1001)
fqz = data.frame(h=fqz$h, f=fqz$f)
fqz$magnitude = abs(fqz$h)
fqz$magnitude_dB = 20 * log10(fqz$magnitude)
fqz$phase = Arg(fqz$h) * 180 / pi

p0 = ggplot(fqz, aes(x = f, y = magnitude_dB)) +
    geom_line(color = "blue", size = 1) +
    geom_vline(xintercept = 50, color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(45, 55), color = "orange", linetype = "dashed") +
    geom_hline(yintercept = -3, color = "darkgreen", linetype = "dashed") +
    labs(title = "Charakterystyka amplitudowa - skala dB",
         x = "Częstotliwość [Hz]", y = "Amplituda [dB]") +
    xlim(30, 70) +
    ylim(-60, 5) +
    theme_minimal()


# Wykres amplitudy w dB
p1 <-ggplot(response, aes(x = frequency, y = magnitude_dB)) +
    geom_line(color = "blue", size = 1) +
    geom_vline(xintercept = 50, color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(45, 55), color = "orange", linetype = "dashed") +
    geom_hline(yintercept = -3, color = "darkgreen", linetype = "dashed") +
    labs(title = "Charakterystyka amplitudowa - skala dB",
         x = "Częstotliwość [Hz]", y = "Amplituda [dB]") +
    xlim(30, 70) +
    ylim(-60, 5) +
    theme_minimal()


# Wykres fazy
p2 <- ggplot(response, aes(x = frequency, y = phase)) +
    geom_line(color = "purple", size = 1) +
    geom_vline(xintercept = 50, color = "red", linetype = "dashed") +
    labs(title = "Charakterystyka fazowa",
         x = "Częstotliwość [Hz]", y = "Faza [stopnie]") +
    xlim(30, 70) +
    theme_minimal()

show(p2)
show(p1)
show(p0)

# Współczynniki filtru
cat("\nWspółczynniki filtru:\n")
cat("Licznik (b):", filter$num, "\n")
cat("Mianownik (a):", filter$den, "\n")