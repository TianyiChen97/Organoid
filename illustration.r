# Define the x and y values for each segment of the function
x <- c(0, 3, 6)
y <- c(0, 2, 6)

# Plot the piecewise linear function with bigger x and y labels
plot(x, y, type = "l", lwd = 4, xlab = "time", ylab = "mirror", cex.lab = 3, mgp = c(1, 1, 0))

# Set the aspect ratio of the plot to be square
asp <- diff(range(x)) / diff(range(y))
plot(x, y, type = "l", lwd = 10, xlab = "time", ylab = "mirror", cex.lab = 5, mgp = c(1, 1, 0), asp = asp)




