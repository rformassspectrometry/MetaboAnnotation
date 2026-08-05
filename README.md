# MetaboAnnotation

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MetaboAnnotation/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MetaboAnnotation/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/MetaboAnnotation/branch/main/graph/badge.svg?token=RLMRUJGOQD)](https://codecov.io/gh/rformassspectrometry/MetaboAnnotation)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/MetaboAnnotation.svg)](https://bioconductor.org/packages/release/bioc/html/MetaboAnnotation.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/MetaboAnnotation.svg)](https://bioconductor.org/packages/stats/bioc/MetaboAnnotation/)
[![build release](http://bioconductor.org/shields/build/release/bioc/MetaboAnnotation.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MetaboAnnotation/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/MetaboAnnotation.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MetaboAnnotation/)

## Overview

High level functions to assist in annotation of (metabolomics) data sets. These
include functions to perform simple tentative annotations based on mass matching
but also functions to consider m/z and retention times for annotation of LC-MS
features given that respective reference values are available. In addition, the
function provides high-level functions to simplify matching of LC-MS/MS spectra
against spectral libraries and objects and functionality to represent and manage
such matched data.

For more information see the package
[homepage](https://rformassspectrometry.github.io/MetaboAnnotation).

## ⤵️ Installation

The package can be installed with

```r
install.packages("BiocManager")
BiocManager::install("MetaboAnnotation")
```

## 🤝 Contribution

Please help us improving and completing the package! Any type of contribution
welcome :open_hands: - including discussions, suggestions or actual code. Don't
be afraid - we're friendly ☺️! 👉 get involved by opening an [issue](https://github.com/rformassspectrometry/MetaboAnnotation/issues).

Please also check out the [**RforMassSpectrometry Contributions
Guide**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).

### 📜 Code of Conduct

We follow the [**RforMassSpectrometry Code of
Conduct**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct)
to maintain an inclusive and respectful community.

## License

This package is licensed under the **Artistic 2.0** license:
📄 [https://opensource.org/license/Artistic-2.0](https://opensource.org/license/Artistic-2.0)

Documentation (manuals, vignettes) is licensed under **CC BY-NC-SA 4.0**:
📄 [https://creativecommons.org/licenses/by-nc-sa/4.0/](https://creativecommons.org/licenses/by-nc-sa/4.0/)
