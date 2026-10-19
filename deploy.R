# ---------------------------------------------------------------------------
# Wavelet Approach - deployment helper. Author: Peterson Owusu Junior.
# One-off helper: publish the app to shinyapps.io (free tier is enough to start)
#
# 1. Create an account at https://www.shinyapps.io and open
#    Account > Tokens > Show, then paste the setAccountInfo() line below.
# 2. Run this script from R (or `Rscript deploy.R`).
# 3. Embed the resulting URL on your website with an <iframe> (see README.md).
# ---------------------------------------------------------------------------
library(rsconnect)

# rsconnect::setAccountInfo(name = "<your-account>", token = "<TOKEN>", secret = "<SECRET>")

app_dir <- normalizePath(dirname(if (interactive()) rstudioapi::getActiveDocumentContext()$path else sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))

rsconnect::deployApp(
  appDir   = app_dir,
  appName  = "wavelet_approach",
  appFiles = c("app.R", list.files("R", full.names = TRUE), list.files("data", full.names = TRUE),
               list.files("www", full.names = TRUE)),
  forceUpdate = TRUE
)
