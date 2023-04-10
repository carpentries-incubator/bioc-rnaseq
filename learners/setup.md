---
title: Setup
---

- Install R, RStudio and packages (see below).

### R and RStudio

- R and RStudio are separate downloads and installations. R is the
  underlying statistical computing environment, but using R alone is
  no fun. RStudio is a graphical integrated development environment
  (IDE) that makes using R much easier and more interactive. You need
  to install R before you install RStudio. After installing both
  programs, you will need to install some specific R packages within
  RStudio. Follow the instructions below for your operating system,
  and then follow the instructions to install packages.

### You are running Windows

<br>

:::::::::::::::  solution

## If you already have R and RStudio installed

- Open RStudio, and click on "Help" > "Check for updates". If a new version is
  available, quit RStudio, and download the latest version for RStudio.

- To check which version of R you are using, start RStudio and the first thing
  that appears in the console indicates the version of R you are
  running. Alternatively, you can type `sessionInfo()`, which will also display
  which version of R you are running. Go on
  the [CRAN website](https://cran.r-project.org/bin/windows/base/) and check
  whether a more recent version is available. If so, please download and install
  it. You can [check here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#How-do-I-UNinstall-R_003f) for
  more information on how to remove old versions from your system if you wish to do so.

- Follow the steps in the instructions [for everyone](#for-everyone) at the
  bottom of this page.
  

:::::::::::::::::::::::::

:::::::::::::::  solution

## If you don't have R and RStudio installed

- Download R from
  the [CRAN website](https://cran.r-project.org/bin/windows/base/release.htm).

- Run the `.exe` file that was just downloaded

- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)

- Under *All Installers* select **RStudio xxxx.yy.zz-uuu.exe - Windows 10/11** (where x, y, z, and u represent version numbers)

- Double click the file to install it

- Once it's installed, open RStudio to make sure it works and you don't get any
  error messages

- Follow the steps in the instructions [for everyone](#for-everyone) at the
  bottom of this page.
  

:::::::::::::::::::::::::

### You are running macOS

<br>

:::::::::::::::  solution

## If you already have R and RStudio installed

- Open RStudio, and click on "Help" > "Check for updates". If a new version is
  available, quit RStudio, and download the latest version for RStudio.

- To check the version of R you are using, start RStudio and the first thing
  that appears on the terminal indicates the version of R you are running. Alternatively, you can type `sessionInfo()`, which will
  also display which version of R you are running. Go on
  the [CRAN website](https://cran.r-project.org/bin/macosx/) and check
  whether a more recent version is available. If so, please download and install
  it.

- Follow the steps in the instructions [for everyone](#for-everyone) at the
  bottom of this page.
  

:::::::::::::::::::::::::

:::::::::::::::  solution

## If you don't have R and RStudio installed

- Download R from
  the [CRAN website](https://cran.r-project.org/bin/macosx/).

- Select the `.pkg` file for the latest R version

- Double click on the downloaded file to install R

- It is also a good idea to install [XQuartz](https://www.xquartz.org/) (needed
  by some packages)

- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)

- Under *All Installers* select **RStudio xxxx.yy.zz-uuu.dmg - macOS 10.15+** (where x, y, z, and u represent version numbers)

- Double click the file to install RStudio

- Once it's installed, open RStudio to make sure it works and you don't get any
  error messages.

- Follow the steps in the instructions [for everyone](#for-everyone) at the
  bottom of this page.
  

:::::::::::::::::::::::::

### You are running Linux

<br>

:::::::::::::::  solution

## Install R using your package manager and RStudio

- Follow the instructions for your distribution
  from [CRAN](https://cloud.r-project.org/bin/linux), they provide information
  to get the most recent version of R for common distributions. For most
  distributions, you could use your package manager (e.g., for Debian/Ubuntu run
  `sudo apt-get install r-base`, and for Fedora `sudo yum install R`), but we
  don't recommend this approach as the versions provided by this are
  usually out of date. In any case, make sure you have at least R 4.2.0.
- Go to the [RStudio download
  page](https://www.rstudio.com/products/rstudio/download/#download)
- Under *All Installers* select the version that matches your distribution, and
  install it with your preferred method (e.g., with Debian/Ubuntu `sudo dpkg -i rstudio-xxxx.yy.zz-uuu-amd64.deb` at the terminal).
- Once it's installed, open RStudio to make sure it works and you don't get any
  error messages.
- Follow the steps in the [instructions for everyone](#for-everyone)
  

:::::::::::::::::::::::::

### For everyone

After installing R and RStudio, you need to install a couple of
packages that will be used during the workshop. We will also learn
about package installation during the course to explain the following
commands. For now, simply follow the instructions below:

- Start RStudio by double-clicking the icon and then type:

```r
install.packages("BiocManager")
BiocManager::install(c("tidyverse", "SummarizedExperiment",
                       "ExploreModelMatrix", "org.Mm.eg.db",
                       "DESeq2", "vsn", "ComplexHeatmap",
                       "RColorBrewer", "hexbin", "cowplot",
                       "clusterProfiler", "enrichplot", "kableExtra"))
```








