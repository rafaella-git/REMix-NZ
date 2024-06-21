---
title: Installation
lang: en-US
---

(installation_label)=

# Installation

In this section you will learn how to set up your local machine to create and
run REMix models. For the installation, you need:

- a Python installation (>=3.9), for example managed with mamba.
- a recent [GAMS](https://www.gams.com/products/gams/gams-language/)
installation (version 37 or above), which is in your PATH in the environment
variables (more detailed explanation below),
  - Make sure your GAMS license file `gamslice.txt` is placed in the GAMS
  installation directory at `C:\GAMS\VERSION#\` or in `~/Documents/GAMS` for
  newer installations.
  - However, in the beginning, a GAMS
  [demo license](https://www.gams.com/try_gams/) is perfectly fine.

```{note}
The REMix framework is written in GAMS but the input data for REMix is mainly
prepared in Python. Therefore, we will install a REMix-API for Python.
Furthermore, we are working with Mamba to create Python environments.
If you are not familiar with this concept, you may refer to this
[user guide](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html).

Additionally, the workflow provided here assumes that you are familiar with Git.
If you are not, there are many useful references online, e.g.
[here](https://www.w3schools.com/git/git_clone.asp?remote=github).
```

We are hosting REMix open source on GitLab. You can clone it from there using
your local Git installation

```bash
git clone https://gitlab.com/dlr-ve/esy/remix/framework.git
```

or by downloading [the repository](https://gitlab.com/dlr-ve/esy/remix/framework)
manually and unzipping the content to a system folder of your choice.

The following steps are shown with Mamba as example, but should work analogously
with other Python environment managers like ``virtualenv`` or ``conda``.

When wanting to install the `remix.framework` package from PyPI, you can do,
e.g., as follows:

```bash
mamba create -n remix-env python=3.10
mamba activate remix-env
pip install remix.framework
```

```{attention}
Please make sure that the installed Python version is supported by your
installed GAMS version. You can check
[on the GAMS website](https://www.gams.com/latest/docs/API_PY_OVERVIEW.html) if
that is the case.
```

You can also create a REMix environment after having cloned or downloaded the
repository to a system directory of your choice. For that, navigate into the
cloned folder and install with `-e` flag:

```bash
cd framework
pip install -e .
```

In case something went wrong, please have a look at the list of frequent
issues below or open an issue on
[our REMix GitLab](https://gitlab.com/dlr-ve/esy/remix/framework/-/issues).

In order to successfully run the tutorials some additional packages are needed.
They can be installed as follows

```bash
pip install remix.framework[tutorials]
```

from PyPI or, after cloning the repository, you can also do

```bash
pip install -e .[tutorials]
```

in your environment and framework directory similar to above.

## Use of GAMS/Python API depending on installed GAMS version

Depending on which version of GAMS you have installed on your machine, different
steps are needed to be executed in your REMix environment.

- GAMS 37-41: no further action needed
- GAMS 42-44: `pip install gams[transfer] --find-links [PATH TO GAMS]\api\python\bdist`
- since GAMS 45: `pip install gamsapi[transfer]==<version_number>.y.z`

## Create and run a REMix model

To create and run REMix models after the installation, you will need to
set up project folders for every model. We provide some tutorials in the
{ref}`tutorial's section <tutorials_basic_label>`, which can be used to gain
first experience working with the REMix optimization framework. Furthermore,
there are some larger example projects available on the
[project's GitLab repository](https://gitlab.com/dlr-ve/esy/remix/framework/-/tree/dev/tutorials).
You can manually download them or clone the repository (cf. above).

To run any of the examples, you can navigate into the project's directory and
use the `remix run` command.

```bash
cd ../project
remix run --datadir=./data/
```

If you want to start your own project, the organization of a project on your
machine should look like this:

```bash
project/
├── data/
|   ├── scen1/
|   |   ├── subscen1/
|   |   |   └── ...
|   |   └── *.dat
|   ├── scen2/
|   └── *.dat
└── ...
```

To learn more about the structure and inheritance logic please check out the
{ref}`inheritance section <explanations_inheritance_label>`.

## Update REMix to the latest version

In your environment run

```bash
pip install remix.framework --upgrade
```

We list notable changes of the API on the changelog.

## Use the developer version of REMix

Checkout the desired branch of REMix and pull the latest changes.
Branches available are:

- main: Latest stable release of REMix. API-breaking changes will be announced.
- dev: Development branch, which is mostly stable. API-breaking changes can be
  included into this branch without any prior notice.
- all others: Specific changes to implement single features or fix bugs.

These branches are typically under active development and are not ensured (see
below to get the list of branches available).

- specific version: You can check out a specific version of REMix, e.g. `0.10.0`.

Check out any branch by using

```bash
git checkout name-of-the-branch-or-version
git pull -r
```

and replacing `name-of-the-branch-or-version` with, for example, `dev`.

## Adding GAMS as environment variable

When running the steps above, you might be running into problems.
One could be that your GAMS installation did not execute correctly.
For that to work, you might have to set system environment variables yourself.
Search in the Windows menu for `Edit environment variables for your account` or
`Umgebungsvariablen für dieses Konto bearbeiten` depending on your language
settings.
In the upper part, create a new user variable called `GAMS_DIR` and give it the
value `C:\GAMS\VERSION#;C:\GAMS\VERSION#\gbin` replacing `VERSION#` with your
installed version of GAMS.

Make sure that `C:\GAMS\VERSION#` is also added to the environment variable
`PATH`.
Usually it is added during installation of GAMS, but sometimes it is not.
You can do that, but just adding the `GAMS_DIR` variable created before.
In Windows, you can do that by adding a line to `Path` containing `%GAMS_DIR%`
(which is referring to the variable between the percent signs).

## Check if GAMS installation is working

To check if GAMS is correctly installed on your machine, you can open a Git
Bash and enter `gams`.
This should prompt something similar to the following links of code:

``` bash
--- Job ? Start 03/24/23 12:16:07 41.1.0 1682d454 WEX-WEI x86 64bit/MS Windows
***
*** GAMS Base Module 41.1.0 1682d454 Oct 28, 2022          WEI x86 64bit/MS Window
***
*** GAMS Development Corporation
*** 2751 Prosperity Ave, Suite 210
*** Fairfax, VA 22031, USA
*** +1 202-342-0180, +1 202-342-0181 fax
*** support@gams.com, www.gams.com
***
*** GAMS Release     : 41.1.0 1682d454 WEX-WEI x86 64bit/MS Windows
*** Release Date     : Oct 28, 2022
*** To use this release, you must have a valid license file for
*** this platform with maintenance expiration date later than
*** Oct 28, 2022
*** System Directory : C:\Program Files\GAMS\41\
***
*** License          : C:\ProgramData\GAMS\gamslice.txt
*** OSMSES Course License                          S230324|0002AB-GEN
*** DLR, Institut fuer Vernetzte Energiesysteme | Energiesystemanalyse
*** DCE2742 01CP
*** License Admin: Name Surname, Email
***
*** Licensed platform                             : Generic platforms
*** The installed license is valid.
*** Evaluation expiration date (GAMS base module) : Apr 23, 2023
*** Note: For solvers, other expiration dates may apply.
*** Run gamslib model licememo for more details.
*** Status: Normal completion
--- Job ? Stop 03/24/23 12:16:07 elapsed 0:00:00.003
```
