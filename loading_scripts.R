library(testthat)
library(devtools)

use_import_from('expm',"%^%")
use_import_from('dplyr',"%>%")
use_import_from('label.switching',"ecr")
use_import_from('stats',"kmeans")
use_import_from('stats',"hclust")
use_import_from('stats',"dist")
use_import_from('stats',"cutree")
use_import_from('truncnorm',"rtruncnorm")
use_import_from('mclust',"Mclust")
use_import_from('mclust',"mclustBIC")

use_import_from("stats", "glm")
use_import_from("stats", "lm")
use_import_from("stats","model.matrix")
use_import_from("stats","predict")
use_import_from("stats", "quasibinomial")
use_import_from("stats", "rbinom")
use_import_from("stats", "rnorm")
use_import_from("stats", "runif")
use_import_from("stats", "sd")

use_test("network_features")
use_test("SBM_estimators")

use_r("assign_treatments")
use_r("seeding")
use_r("optimal_design")


document()


