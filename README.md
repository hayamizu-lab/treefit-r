# Treefit for R

## Install

### Debian GNU/Linux, Ubuntu

```bash
sudo -H apt install -V -y r-base libcurl4-openssl-dev
sudo -H Rscript -e 'install.packages(c("devtools"))'
sudo -H Rscript -e 'devtools::install_github("hayamizu-lab/treefit-r")'
```

### macOS

```bash
brew cask install r
echo 'options(repos="https://cloud.r-project.org")' >> ~/.Rprofile
Rscript -e 'install.packages(c("devtools"))'
sudo -H Rscript -e 'devtools::install_github("hayamizu-lab/treefit-r")'
```

## Usage

```R
library(treefit)
tree <- generate_2d_n_wands_expression(500, 5, 0.1)
plot(tree)
estimated <- estimate(list(expression=tree))
plot_estimated(estimated)
```
