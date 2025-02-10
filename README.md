# Treefit for R - The first software for quantitative trajectory inference

This is an implementation of
[**Treefit**](https://hayamizu-lab.github.io/treefit/) in R.

**Treefit** is a novel data analysis toolkit that helps you perform
two types of analysis: 1) checking the goodness-of-fit of tree models
to your single-cell gene expression data; and 2) deciding which tree
best fits your data. Treefit for R can be used in conjunction with
other popular software packages, such as
[Seurat](https://satijalab.org/seurat/) and
[dynverse](https://dynverse.org/).

We'll implement [Monocle
3](https://cole-trapnell-lab.github.io/monocle3/) integration soon.

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

The main functions are `treefit::treefit()` and `plot()`:

```R
fit <- treefit::treefit(YOUR_SINGLE_CELL_GENE_EXPRESSION_DATA)
plot(fit)
```

See `vignette("treefit")` for details.

## How to release

Bump `Version:` in `DESCRIPTION`.

Add a release note to `NEWS.md`. You need to run `Rscript -e
"devtools::spell_check()` after it. You may need to update
`inst/WORDLIST` too.

Check on local:

```bash
Rscript -e "roxygen2::roxygenise()"
Rscript -e "rcmdcheck::rcmdcheck(args=c('--as-cran'), error_on='warning', check_dir='check')"
```

Check on win-builder.r-project.org:

```bash
rake check:win_builder
```

Update `cran-comments.md`.

Submit to CRAN:

```bash
rake cran:submit
```

Release after the submission is accepted:

```bash
rake release
```
