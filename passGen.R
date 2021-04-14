
paste(sample(c(0:9, letters, LETTERS), 24, replace = T), collapse = "")



len = 16e3
len = 10e3
len = 6e3
text = sample(letters, len, replace = T)
text[sample(1:len, len * .15)] = ' '
paste0(text, collapse = '')


