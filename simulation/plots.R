plot_risk <- function(){

    # box plot risk by protection level
p1 <- ggplot(individuals, aes(x = as.factor(individual_protected), y = individual_risk)) + # create box plot of risk by protection level
    geom_boxplot() +
    labs(x = "Protection Level", y = "Risk") 

# box plot risk by size class
p2 <- ggplot(individuals, aes(x = as.factor(size_individual), y = individual_risk)) + # create box plot of risk by size class
    geom_boxplot() +
    labs(x = "Size Class", y = "Risk") 

# box plot risk by foraging guild
p3 <- ggplot(individuals, aes(x = as.factor(foraging_individual), y = individual_risk)) + # create box plot of risk by foraging guild
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Risk") 

# historgram of risk

p4 <- ggplot(individuals, aes(x = individual_risk)) + # create histogram of risk
    geom_density() +
    labs(x = "Risk", y = "Count") 

(p1 + p2 + p3 )/ p4# combine plots

}

# plot behaviours

plot_beh <- function(){
    # box plot behaviours by protection level
p1 <- ggplot(response, aes(x = as.factor(individual_protected), y = foraging)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Foraging")

p2 <- ggplot(response, aes(x = as.factor(individual_protected), y = vigilance)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Vigilance")

p3 <- ggplot(response, aes(x = as.factor(individual_protected), y = movement)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Movement")

p4 <- ggplot(response, aes(x = as.factor(individual_protected), y = bites)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Bites")

# box plot behaviours by size class
p5 <- ggplot(response, aes(x = as.factor(size_individual), y = foraging)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Foraging")

p6 <- ggplot(response, aes(x = as.factor(size_individual), y = vigilance)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Vigilance")

p7 <- ggplot(response, aes(x = as.factor(size_individual), y = movement)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Movement")

p8 <- ggplot(response, aes(x = as.factor(size_individual), y = bites)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Bites")

# box plot behaviours by foraging guild
p9 <- ggplot(response, aes(x = as.factor(foraging_individual), y = foraging)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Foraging")

p10 <- ggplot(response, aes(x = as.factor(foraging_individual), y = vigilance)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Vigilance")

p11 <- ggplot(response, aes(x = as.factor(foraging_individual), y = movement)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Movement")

p12 <- ggplot(response, aes(x = as.factor(foraging_individual), y = bites)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Bites")

# Combine plots grouped by predictor
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8  + p9 + p10 + p11 + p12 + plot_layout(ncol = 4, nrow = 3, byrow = TRUE) 

}

beh_dis <- function(){
    p1 <- ggplot(response, aes(x = movement)) +
    geom_density() +
    labs(x = "Movement", y = "Count") 

p2 <- ggplot(response, aes(x = foraging)) +
    geom_density() +
    labs(x = "Foraging", y = "Count")

p3 <- ggplot(response, aes(x = vigilance)) +
    geom_density() +
    labs(x = "Vigilance", y = "Count")

p4 <- ggplot(response, aes(x = bites)) +
    geom_density() +
    labs(x = "Bites", y = "Count")

p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow = 2, byrow = TRUE) # combine plots
}