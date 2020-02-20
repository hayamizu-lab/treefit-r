# Treefit for R - Quantitative analysis software for trajectory inference

This is an implementation of [**Treefit**](https://hayamizu-lab.github.io/treefit/) in R.

**Treefit** is a novel data analysis toolkit that helps you perform two types of analysis: 1) checking the goodness-of-fit of tree models to your single-cell gene expression data; and 2) deciding which tree best fits your data. Treefit for R can be used in conjunction with other popular software packages, such as Monocle, Seurat, and dynverse.

## Install

### Debian GNU/Linux and Ubuntu

```bash
sudo -H apt install -V -y r-base libcurl4-openssl-dev
sudo -H Rscript -e 'install.packages(c("treefit"))'
```

### macOS

```bash
brew cask install r
echo 'options(repos="https://cloud.r-project.org")' >> ~/.Rprofile
Rscript -e 'install.packages(c("treefit"))'
```

### Windows

```r
install.packages(c("treefit"))
```

## Usage

The main functions are `treefit::estimate()` and `treefit::plot_estimated()`:

```R
estimated <- treefit::estimate(YOUR_SINGLE_CELL_GENE_EXPRESSION_DATA)
treefit::plot_estimated(estimated)
```

See the following vignettes for more details:

  * `vignette("treefit")`
  * `vignette("seurat-integration")`
