# ---------------------------------------------------------------------------
# Wavelet Approach - deployment helper. Author: Peterson Owusu Junior.
# One-off helper: publish the app to shinyapps.io (free tier is enough to start)
#
# Account "pete-econ" is registered on this machine (rsconnect::accounts()).
# On a new machine: shinyapps.io > Account > Tokens > Show, then run
#   rsconnect::setAccountInfo(name = "pete-econ", token = "<TOKEN>", secret = "<SECRET>")
# Run this script from R (source("deploy.R")) or `Rscript deploy.R`.
# NOTE: the app now lives on Posit Connect Cloud (auto-published from GitHub);
# the shinyapps.io copy was retired on 2026-10-20. This script only matters if
# you want to bring the shinyapps.io copy back.
# ---------------------------------------------------------------------------
library(rsconnect)

app_dir <- tryCatch(
  normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))),
  error = function(e) NA)
if (is.na(app_dir) || !file.exists(file.path(app_dir, "app.R"))) app_dir <- "C:/Users/peteo/wavelet/wavelet_app"

rsconnect::deployApp(
  appDir   = app_dir,
  appName  = "wavelet_approach",
  account  = "pete-econ",
  appFiles = c("app.R", list.files("R", full.names = TRUE), list.files("data", full.names = TRUE),
               list.files("www", full.names = TRUE)),
  forceUpdate = TRUE
)
