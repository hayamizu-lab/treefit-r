## Resubmission

This is a resubmission.

Previous submission was rejected by auto pre-tests on win-builder:

> package treefit_1.0.1.tar.gz does not pass the incoming checks automatically, please see the following pre-tests:
> Windows: <https://win-builder.r-project.org/incoming_pretest/treefit_1.0.1_20210319_233920/Windows/00check.log>
> Status: 1 WARNING
> Debian: <https://win-builder.r-project.org/incoming_pretest/treefit_1.0.1_20210319_233920/Debian/00check.log>
> Status: OK

The warning is caused by download error in a vignette. The vignette
downloads 6MB data from the Internet. We've fixed the warning by just
ignoring the download error because the vignette is an optional
vignette that depends on a package in Suggests, Seurat.

## Test environments

* local: R 4.0.4 x86_64-pc-linux-gnu (64-bit)
* GitHub Actions:
  * R 4.0 on Ubuntu 20.4
  * R 4.0 on macOS Catalina 10.15
  * R 4.0 on Windows Server 2019
  * R devel on Windows Server 2019
* win-builder:
  * R 4.0 on x86_64-w64-mingw32
  * R devel on x86_64-w64-mingw32

## R CMD check results

Status: OK

## revdepcheck results

There are currently no downstream dependencies for this package.
