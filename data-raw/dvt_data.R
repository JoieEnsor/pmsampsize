## Data preparation script for dvt_data
## Source: binary_data_example.csv (simulated DVT dataset, n=9805)
## Run this script to regenerate data/dvt_data.rda

dvt_data <- read.csv("../binary_data_example.csv")
dvt_data <- dvt_data[, c("male", "OCuse", "malignancy", "surgery",
                           "ab_leg_trauma", "vein_distension", "calf_diff",
                           "ab_ddimer", "pregnancy", "immobile", "clotting",
                           "dvt")]
save(dvt_data, file = "data/dvt_data.rda", compress = "bzip2")
