library(neuralnet)
data_xor <- data.frame(
    input1 = c(0,0,1,1),
    input2 = c(0,1,0,1),
    output = c(0,1,1,0)
)

nn1 <- neuralnet(output ~ input1 + input2,
                 data = data_xor,
                 hidden = c(2),
                 linear.output = FALSE)

nn2 <- neuralnet(output ~ input1 + input2,
                 data = data_xor,
                 hidden = c(2),
                 linear.output = TRUE)

pred1 <- compute(nn1, data_xor[,1:2])$net.result
pred2 <- compute(nn2, data_xor[,1:2])$net.result

results <- data.frame(
    Expected = data_xor$output,
    Nonlinear = round(pred1, 3),
    Linear = round(pred2, 3)
)
print(results)