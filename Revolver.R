library(devtools)

devtools::install_github("caravagn/revolver")

library(revolver)
library(dplyr)

revolver_data <- read.csv("/media/kovac/Resources1/MartinD/revolver_input.tsv",sep = "\t")
view(revolver_data)
revolver_data["cluster"] <- as.character(revolver_data["cluster"])
revolver_data["is.driver"] <- as.logical(unlist(revolver_data["is.driver"]))
revolver_data["is.clonal"] <- as.logical(unlist(revolver_data["is.clonal"]))
data_cohort <- revolver_cohort(revolver_data)
data_cohort_copy <- data_cohort

data_cohort$variantIDs <- unique(data_cohort$variantIDs)
View(data_cohort$variantIDs)


# Remove duplicated values
data_cohort$patients


revolver::revolver_check_cohort(data_cohort,stopOnError = F)



revolver_fit <- revolver::revolver_fit(trees,)

revolver_burden_fit <- revolver::revolver_fit(trees_burden)

trees <- revolver::compute_clone_trees(data_cohort)
trees_burden <- revolver::compute_mutation_trees(data_cohort)
revolver::revolver_check_cohort(revolver_fit)
revolver::plot_trajectories_per_cluster(data_cohort)
View(trees)
#plot_patient_trees(trees,"P")
revolver::plot_patient_mutation_burden()

