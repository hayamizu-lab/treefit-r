## Resubmission

This is a resubmission.

> The LICENSE file is only needed if you have additional restrictions to
> the GPL-3 which you have not? In that case omit the file and its
> reference in the DESCRIPTION file.
> "GPL-3 | file LICENSE" --> "GPL-3"
> If you do have additional restrictions that I failed to see, please
> change the reference to "GPL-3 + file LICENSE" and change the file
> LICENSE to *only* contain the additional restrictions.

We don't have additional restrictions. So we've changed "GPL (>= 3) |
file LICENSE" to "GPL (>= 3)".

> message() is easily suppressed. You don't need to ask the user for
> permission to use it.

We've removed the "verbose" argument to control whether calling
message() or not. And we've also removed needless message()s.

> Please add \value to all .Rd files that are not data files and explain
> the functions results in the documentation.
> f.i.: plot.treefit.Rd
> If a function does not return a value, please document that too,
> e.g. \value{None}.

We've added \value to all .Rd files. We don't have any functions that
don't return a value. So we've described a return value for all
functions.

### Other changes since the first submission

We have the following additional changes since the first submission:

  * Improved documents. (No API change.)

  * Fixed a bug in treefit::treefit() function. (No API change.)

    * Added a test for the case.

  * Added a new "n_perturbations" option to treefit::treefit() function.
    It's documented. This changes API.

  * Improved internal computation in treefit::treefit() function.
    (No API change.)

  * Improved visualization by treefit::plot.treefit() function.
    (No API change.)

## Test environments

* local: x86_64-pc-linux-gnu-3.6.2
* GitHub Actions:
  * R 3.5 on Ubuntu 18.04
  * R 3.6 on Ubuntu 18.04
  * R 3.6 on macOS Catalina 10.15
  * R 3.5 on Windows Server 2019
  * R 3.6 on Windows Server 2019

## R CMD check results

Status: 1 NOTE

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Kouhei Sutou <kou@clear-code.com>’

    New submission

This NOTE is acceptable because this is a new package.

## revdepcheck results

There are currently no downstream dependencies for this package.
