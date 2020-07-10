## Gaussian blur
data <- read.csv("/Users/chaokuan-hao/Documents/BIO_IT_Station/Proteomics/alldata.ctrl.csv", header=T, sep=",")
colnames(data)

H_4_row <- data[grep("Percentage.H.", colnames(data),value=T)][4,]
H_4_row

plot(c(1:27), H_4_row)

myMatrix <- as.matrix(H_4_row)

# Create a sequence of numbers between -10 and 10 incrementing by 0.1.
x <- seq(1, 11, by = 1)
# Choose the mean as 2.5 and standard deviation as 0.5.
y <- dnorm(x, mean = 6, sd = 0.5)
# plot(x,y)

kernel <- as.matrix(y)
convolved_result <- convolve(myMatrix, kernel, type="open")
plot(c(1:27), convolved_result[6:32])

## PrInCE


