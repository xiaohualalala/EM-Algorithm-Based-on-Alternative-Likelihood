library(ggplot2)
library(caret)

mixNB_density <- function(x, sizes, probs, weights) {
    # Initialize the density vector
    density <- numeric(length(x))
    # Compute the mixed density
    for (i in 1:length(x)) {
        for (j in 1:length(weights)) {
            # Add the weighted density of each Negative Binomial component
            density[i] <- density[i] + weights[j] * dnbinom(x[i], size = sizes[j], prob = probs[j])
        }
    }
    return(density)
}

plotHistDensity <- function(observations, estimation, accuracy) {
    eSize <- estimation$eSize
    eProb <- estimation$eProb
    eWeight <- estimation$eWeight

    # Generate and plot the density estimation for the observations
    range <- seq(min(observations), max(observations))
    density <- mixNB_density(range, eSize, eProb, eWeight)
    data <- data.frame(Value = observations)
    density_data <- data.frame(x = range, density = density)
    # Plotting using ggplot2
    p <- ggplot(data, aes(x = Value)) +
        geom_histogram(aes(y = ..density..), binwidth = diff(range(observations)) / 55, fill = "skyblue", color = "white") +
        geom_line(data = density_data, aes(x = x, y = density), color = "red", size = 1) +
        theme_bw() +
        labs(x = "Value", y = "Density") + 
        annotate("text", x = Inf, y = Inf, label = paste("Accuracy:", accuracy*100, "%"), hjust = 1.1, vjust = 2, size = 4, color = "blue")

    # Display the plot
    print(p)
}

