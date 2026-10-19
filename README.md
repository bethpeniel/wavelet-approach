# Wavelet Approach – Bivariate & Multiple Wavelet Analysis (Shiny app)

**Author:** Peterson Owusu Junior, University of Cape Coast, Ghana
— <peterson.owusu@ucc.edu.gh> · <powusujunior@gmail.com> · <https://positivecarryacademy.com/>
([Google Scholar](https://scholar.google.com/citations?user=r6ZueU4AAAAJ&hl=en) ·
[Scopus](https://www.scopus.com/authid/detail.uri?authorId=58833238800) ·
[LinkedIn](https://www.linkedin.com/in/peterson-owusu-junior-phd-80858334/) ·
[ResearchGate](https://www.researchgate.net/profile/Peterson-Owusu-Junior) ·
[IDEAS/RePEc](https://ideas.repec.org/e/pow25.html))

A browser-based app, modelled on David Gabauer's *Connectedness Approach* app, for
studying co-movement of financial and economic time series in the time–frequency domain.

Suggested citation: Owusu Junior, P. (2026). *Wavelet Approach: A Shiny application for
bivariate and multiple wavelet analysis*. R Shiny application.

| Approach | Engine | Outputs |
|---|---|---|
| **Bivariate wavelets** (continuous Morlet) | `WaveletComp` | wavelet power, cross-wavelet power, wavelet coherence with phase arrows, average power/coherence, phase difference (band + image), reconstruction, tables |
| **Multiple wavelets** (MODWT) | `waveslim`, `wavemulcor` (after Fernández-Macho 2012 / Polanco-Martínez's `W2CWM2C`) | MODWT decomposition, pairwise wavelet correlation, wavelet cross-correlation, wavelet multiple correlation (WMC), wavelet multiple cross-correlation (WMCC), wavelet local multiple correlation (WLMC), tables |

Common tabs: data table, raw data, transformed data, summary statistics, and a
**Download results** button that bundles CSV tables and PNG figures into a zip.

## Folder layout

```
wavelet_app/
├── app.R                    UI + server
├── R/
│   ├── data_utils.R         file reading, date parsing, transformations, summary stats
│   ├── bivariate_wavelets.R WaveletComp wrappers (coherence, phase, reconstruction)
│   └── multiple_wavelets.R  MODWT, WC, WCC, WMC, WMCC, WLMC + ggplot2 plots
├── data/                    replication datasets (AFRICA_week, gold_cryptos9, REIT, Ghana cocoa/gold/FX)
├── www/custom.css
├── deploy.R                 helper for shinyapps.io
└── README.md
```

## Run locally

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "WaveletComp", "waveslim",
                   "wavemulcor", "shinycssloaders", "readxl", "zip"))
shiny::runApp("C:/Users/peteo/wavelet/wavelet_app")
```

or from a terminal (Rscript is not on the PATH on this machine, so use the full path):

```powershell
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" -e "shiny::runApp('C:/Users/peteo/wavelet/wavelet_app', port = 4321, launch.browser = TRUE)"
```

Quickest option: double-click `start_app.bat`. It starts the server on port 4321 (if it is not
already running) and opens `open_app.html`, a plain HTML page that embeds the app in a
full-window iframe. If the server is not running, that page shows instructions instead of a
blank frame. `open_app.html` is also the template for embedding the app in a website: change
the `APP_URL` constant at the top of the file to the published address.

## Input data (your own files)

The app is built around uploading new datasets; the bundled example files only illustrate the format.

* Accepted: csv, txt/tsv (comma, semicolon, tab or pipe separated; decimal point or comma;
  thousands separators tolerated) and xlsx/xls (first sheet).
* One column per series, one row per observation, column names in the first row
  (untick *First row contains column names* otherwise).
* A date column is recommended. Pick it in *Date column* (defaults to the first column) or choose
  *No date column* to index observations 1..N. Formats such as `YYYY-MM-DD`, `DD/MM/YYYY`,
  `MM/DD/YYYY`, `DD.MM.YYYY`, `YYYYMMDD`, `Jan-15`, Excel dates and Excel serial numbers are
  auto-detected; if detection is wrong, set *Date format* explicitly.
* Rows are sorted chronologically, rows with missing values are dropped, and non-numeric columns
  are ignored. The app reports what it changed above the result tabs.
* Any number of series is supported. For more than 7 series the pairwise correlation tab switches
  to a per-level correlation matrix layout.
* Upload limit is 30 MB (adjust `shiny.maxRequestSize` in `app.R` if needed).

## Settings worth knowing

* **Transformation**: use *nth log returns* (or percentage changes) for price levels; the
  Gold/crypto and REIT replication sets are already returns.
* **Bivariate**: `dj` controls frequency resolution (1/20 is quick; 1/250 is very slow).
  Monte-Carlo significance with 100 simulations takes a few seconds for ~1,000 observations.
  The *Phase difference* tab lets you pick the period band to average over.
* **Multiple**: the maximum number of levels `J` depends on the sample size and filter
  (boundary coefficients are removed with a brick wall; the app clamps `J` and tells you).
  `Max lag` sets the lead/lag range for cross-correlations. The local window `M` for WLMC
  defaults to N/8.
* Level `dj` covers periods `2^j` to `2^(j+1)` observations (d1 = 2–4, d2 = 4–8, …).

## Source code on GitHub

The app lives at <https://github.com/bethpeniel/wavelet-approach> (private). This folder is a
git clone linked to it, so after editing any file either double-click `push_update.bat` or run:

```powershell
git add -A
git commit -m "describe the change"
git push
```

To make the repository public later: `gh repo edit bethpeniel/wavelet-approach --visibility public`.

## Publish and embed on your website

1. Deploy to shinyapps.io (free tier works):

   ```r
   install.packages("rsconnect")
   rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
   source("C:/Users/peteo/wavelet/wavelet_app/deploy.R")
   ```

   The app is then available at `https://<account>.shinyapps.io/wavelet_approach/`.

2. Embed it in any web page:

   ```html
   <iframe src="https://<account>.shinyapps.io/wavelet_approach/"
           style="width:100%; height:100vh; border:0;" loading="lazy"
           title="Wavelet Approach"></iframe>
   ```

   For WordPress, Wix, Squarespace or similar, paste that snippet into an HTML / embed block.

Alternatives: Posit Connect Cloud, or a self-hosted Shiny Server / ShinyProxy on a VPS.
(`shinylive` does not currently work for this app because `WaveletComp` and `waveslim`
need compiled code that is not available in the WebAssembly R build.)

## References

### Applications by Peterson Owusu Junior

* Owusu Junior, P., Adam, A. M., & Tweneboah, G. (2017). Co-movement of real exchange rates in the West African Monetary Zone. *Cogent Economics & Finance*, 5(1), 1351807. https://doi.org/10.1080/23322039.2017.1351807
* Owusu Junior, P., Boafo, B. K., Awuye, B. K., Bonsu, K., & Obeng-Tawiah, H. (2018). Co-movement of stock exchange indices and exchange rates in Ghana: A wavelet coherence analysis. *Cogent Business & Management*, 5(1), 1481559. https://doi.org/10.1080/23311975.2018.1481559
* Owusu Junior, P., Tweneboah, G., & Adam, A. M. (2019). Interdependence of major exchange rates in Ghana: A wavelet coherence analysis. *Journal of African Business*, 20(3), 407–430. https://doi.org/10.1080/15228916.2019.1583973
* Tweneboah, G., Owusu Junior, P., & Oseifuah, E. K. (2019). Integration of major African stock markets: Evidence from multi-scale wavelets correlation. *Academy of Accounting and Financial Studies Journal*, 23(6).
* Owusu Junior, P., & Tweneboah, G. (2020). Are there asymmetric linkages between African stocks and exchange rates? *Research in International Business and Finance*, 54, 101245. https://doi.org/10.1016/j.ribaf.2020.101245
* Owusu Junior, P., Alagidede, I., & Tweneboah, G. (2020). Shape-shift contagion in emerging markets equities: Evidence from frequency- and time-domain analysis. *Economics and Business Letters*, 9(3), 146–156. https://doi.org/10.17811/ebl.9.3.2020.146-156
* Owusu Junior, P., Adam, A. M., Asafo-Adjei, E., Boateng, E., Hamidu, Z., & Awotwe, E. (2021). Time-frequency domain analysis of investor fear and expectations in stock markets of BRIC economies. *Heliyon*, 7(10), e08211. https://doi.org/10.1016/j.heliyon.2021.e08211
* Asafo-Adjei, E., Boateng, E., Isshaq, Z., Idun, A. A.-A., Owusu Junior, P., & Adam, A. M. (2021). Financial sector and economic growth amid external uncertainty shocks: Insights into emerging economies. *PLOS ONE*, 16(11), e0259303. https://doi.org/10.1371/journal.pone.0259303
* Agyei, S. K., Adam, A. M., Bossman, A., Asiamah, O., Owusu Junior, P., Asafo-Adjei, R., & Asafo-Adjei, E. (2022). Does volatility in cryptocurrencies drive the interconnectedness between the cryptocurrencies market? Insights from wavelets. *Cogent Economics & Finance*, 10(1), 2061682. https://doi.org/10.1080/23322039.2022.2061682
* Agyei, S. K., Owusu Junior, P., Bossman, A., Asafo-Adjei, E., Asiamah, O., & Adam, A. M. (2022). Spillovers and contagion between BRIC and G7 markets: New evidence from time-frequency analysis. *PLOS ONE*, 17(7), e0271088. https://doi.org/10.1371/journal.pone.0271088
* Barson, Z., Owusu Junior, P., & Adam, A. M. (2023). Comovement between commodity returns in Ghana: The role of exchange rates. *Journal of Economic Structures*, 12, 17. https://doi.org/10.1186/s40008-023-00312-z
* Owusu Junior, P., Tetteh, J. E., Nkrumah-Boadu, B., & Adjei, A. N. K. (2024). Comovement of African stock markets: Any influence from the COVID-19 pandemic? *Heliyon*, 10(9), e29409. https://doi.org/10.1016/j.heliyon.2024.e29409
* Yeboah, S. D., Owusu Junior, P., Idun, A. A.-A., Akorsu, P. K., & Fumey, M. P. (2025). Dynamic interdependence of major currencies and the US dollar: A wavelet coherence approach. *Future Business Journal*, 11(1). https://doi.org/10.1186/s43093-025-00636-1
* Tetteh, J. E., & Owusu Junior, P. (2025). Co-movement between stock markets in advanced economies and Africa in times of uncertainty: A time-frequency domain approach. *PLOS ONE*, 20(11), e0334325. https://doi.org/10.1371/journal.pone.0334325

### Software

Peterson Owusu Junior is acknowledged in the W2CWM2C package documentation (`?W2CWM2C`) for reporting bugs that were fixed in the package.

* Roesch, A., & Schmidbauer, H. (2025). *WaveletComp: Computational wavelet analysis*. R package version 1.2. https://doi.org/10.32614/CRAN.package.WaveletComp
* Polanco-Martínez, J. M. (2021). *W2CWM2C: A graphical tool for wavelet (cross) correlation and wavelet multiple (cross) correlation analysis*. R package version 2.2. https://doi.org/10.32614/CRAN.package.W2CWM2C
* Polanco-Martínez, J., & Fernández-Macho, J. (2014). The package 'W2CWM2C': Description, features and applications. *Computing in Science & Engineering*, 16(6), 68–78. https://doi.org/10.1109/MCSE.2014.96
* Fernández-Macho, J. (2021). *wavemulcor: Wavelet routines for global and local multiple regression and correlation*. R package version 3.1.2. https://doi.org/10.32614/CRAN.package.wavemulcor
* Whitcher, B. (2024). *waveslim: Basic wavelet routines for one-, two-, and three-dimensional signal processing*. R package version 1.8.5. https://doi.org/10.32614/CRAN.package.waveslim

### Methodology

* Fernández-Macho, J. (2012). Wavelet multiple correlation and cross-correlation: A multiscale analysis of Eurozone stock markets. *Physica A*, 391(4), 1097–1104. https://doi.org/10.1016/j.physa.2011.11.002
* Fernández-Macho, J. (2018). Time-localized wavelet multiple regression and correlation. *Physica A*, 492, 1226–1238.
* Gençay, R., Selçuk, F., & Whitcher, B. (2001). *An introduction to wavelets and other filtering methods in finance and economics*. Academic Press.
* Whitcher, B., Guttorp, P., & Percival, D. B. (2000). Wavelet analysis of covariance with application to atmospheric time series. *Journal of Geophysical Research: Atmospheres*, 105(D11), 14941–14962.
* Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis. *Bulletin of the American Meteorological Society*, 79, 61–78.
* Ramsey, J. B. (2002). Wavelets in economics and finance: Past and future. *Studies in Nonlinear Dynamics & Econometrics*, 6(3).
