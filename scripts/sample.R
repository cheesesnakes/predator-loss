# creating a subsample of the quadrat videos

## get data on deployments

sites <- read.csv("../data/chapter-3_site_metadata.csv")

## make a data frame with deployment-id

sample <- data.frame(deployment_id = rep(sites$deployment_id, 50))

## add plot id 1 - 5 for each deployment

sample$plot_id <- rep(1:5, 13)

## for each deployment make 10 samples between 600 and 3600 seconds

sample$sample_s <- sample(600:3600, 650, replace = TRUE)

## add a column with the time in minutes

sample$sample_m <- sample$sample_s / 60

## add sample ID 1-10 for each deployment

sample$sample_id <- rep(1:10, 13)

sample

## save the data

write.csv(sample, "chapter-3_samples.csv", row.names = FALSE)
