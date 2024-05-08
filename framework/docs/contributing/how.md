---
title: How to contribute?
lang: en-US
---

(contributing_how)=

# How to contribute?

You want to contribute to the development of REMix?
Welcome!
Below you find some information on how you can do that most easily.

## Workflow for active development

We use `git` as software versioning tool to maintain and develop
REMix.
You can make changes to the source code on your local copy of the
software and create a merge request to include your changes in the
official software.
Information on
[merge requests](https://docs.gitlab.com/ee/user/project/merge_requests)
are provided in the offical GitLab documentation.

### Get the developer version

Before creating a merge request on REMix you have to

1. create a fork of the repository in your own account and
2. delete the `framework` folder containing the REMix installation
   from your machine.
   Note, you do NOT need to reinstall REMix after
   this process.

After this, open a terminal or git bash in the `remix` folder (the
folder where your `framework` installation used to be) and clone your
fork:

```bash
git clone https://gitlab.com/YOUR_GITLAB_USERNAME/framework.git
```

or if you are using ssh:

```bash
git clone git@gitlab.com:YOUR_GITLAB_USERNAME/framework.git
```

This creates a new folder `framework` and downloads the contents from your fork of the source code.
Then install REMix develop and/or documentation requirements:

```bash
cd framework
mamba activate remix-env
pip install -e .[dev]
pip install -e .[docs]
```

In order to stay in sync with the original `remix/framework` repository, we
add a link as remote to your local copy of REMix.
We will call the link `upstream`.
Fetch branches available and after that, you can pull changes
from a specific branch of the REMix repository (branch `dev` in the example
below).
Adjust the link in case you are using `ssh` instead of `https`.

```bash
git remote add upstream https://gitlab.com/dlr-ve/remix/framework.git
git fetch upstream
# checkout your local dev branch and pull updates to from the dev branch of
# the upstream repository (remix/framework) and rebase
git checkout dev
git pull upstream dev --rebase
```

### Create a new branch

We use the `--rebase` command to avoid merge commits for every upstream pull.
If you want to make changes to REMix, you can now follow these steps:

-  Rebase your local `dev` branch from the REMix `dev` branch as described above.
-  Checkout a new branch from your local `dev` branch.
-  In the REMix development, we follow a naming convention for feature branches.
   The name of a should start with
   -  `fix/` if you want to fix something in the code.
   -  `feature/` if you want to add a new feature.
   -  `pages/` if you want to change documentation (this will enable it
      appearing in the public documentation automatically).
-  It is helpful to include the number of the issue the branch/merge request is
   addressing, e.g. `#647`.
-  An example for the creation of a feature branch following this best
   practices could look like this:

   `git checkout -b fix/#647-what-does-your-branch-fix`.

-  Change, add or remove code and/or documentation.
-  Make your changes and commit them.
-  Push your changes to your fork of the repository.
-  Open a new
   [merge request](https://docs.gitlab.com/ee/user/project/merge_requests).
-  Check back for answers of the maintainers.

### Open a merge request

When you start a new merge request, you will be prompted for several
information:

1. Define which branch of your fork should be merged into which branch of the
   framework.
2. Provide the information required by GitLab:

   -  Give a meaningful `Title` to the merge request, which briefly summarizes
      the changes applied.
      If it is still work in progress, you can mark it a draft by adding
      `Draft:` as prefix for the title.
   -  In the `Description` you should provide a short description on what the
      new code does, what additional input files are needed or whih new symbols
      are added to the gdx output file.
      This will help the reviewers to better understand your changes.
   -  For `Assignee` put at least one of the REMix release managers.
   -  For `Reviewers` put someone who also works with REMix and is familiar with
      either your or similar topics (e.g. transfer networks).
   -  `Milestones` and `Labels` are currently not required.

3. After submitting the merge request GitLab will automatically try to
   temporarily merge the code and run the testcases. You should get
   either the information that the code can be merged or that there are
   merge conflicts.
4. If your code can be merged and you verified all of your changes,
   remember to change the title of the merge request and remove the
   `Draft:` prefix to indicate you want to proceed.
5. Reviewers and release managers will take a look at your code and may
   comment about specific lines of code or may ask additional questions.
   If you started the merge request, GitLab will automatically keep you
   informed by mail about any changes.
6. As soon as the changes have been approved by both the reviewer and
   the release manager, the release manager will merge the branch to
   main - now your changes have improved the REMix model for the future!

````{note}
**Tips on resolving merge conflicts**

If you encounter merge conflicts, the best way to handle them is to
locally pull your main branch to an up-to-date state, switch to your
development branch and merge the main branch onto your development
branch.

```bash
git checkout main
git pull -r upstream main
git checkout <new_branch_name>
git merge --no-ff main
```

Now you will get a notification in which files the merge conflicts have
occured and have to decide for each conflict which version of the code
from either your branch or the main branch is the correct one. If you
resolved the merge conflict in a file, you need to manually confirm it
using the `git add` command. After you have resolved all merge
conflicts you can commit your code again and push it to remote - your
merge request will automatically update itself and ideally show you that
it can now be merged.

To verify what changes are applied to the code, it is also quite helpful
to check the `"changes"` tab in your merge request to verify you have
only changed files you wanted to change.

**Where to find what code in REMix framework**

After getting to know how to collaborate, we will quickly show where you
can make changes in REMix. The graph below gives you an overview of the
repository structure:

```bash
framework/
├── docs/
├── remix/framework/
└── testing/
```
````

### Create software tests

Lastly, to ensure changes made in the code do not break with existing
code and to verify the implemented features provide results as expected,
automatic software tests are carried out. This helps to make future
development, application and understanding of the code better. If you
make changes to REMix core components or in the pre- and postprocessing,
we encourage you to write testcases and a short documentation about
options, parameters and how to use them properly.

A merge request will automatically run these tests. However, you can run
the tests on your local machine as well prior to opening the merge
request. To do that, in your Python environment with REMix installed
install the `pytest` package and then run

```bash
python -m pytest
```

from the root of the `framework` code.

If you want to add your own tests, you can add them in the testing subdirectory.
In case you have any questions on setting up specific tests, do not hesitate to
reach out to the developers.

## Tips on contribution for the documentation

If you want to contribute to the documentation, you need to install the docs
dependencies first.
You can do so by running

```bash
pip install -e .[docs]
```

in the `remix/framework` directory and in your local REMix environment.

Most of the documentation is written as Markdown code.
For a good guide on the basic syntax of Markdown see
[https://www.markdownguide.org/basic-syntax](https://www.markdownguide.org/basic-syntax).

### General overview

To get a rough idea on which type of documentation to place where, the following
table might be helpful:

| Section on GitLab pages (docu type) | Purpose                                          | Target group                                          | What / How to document                                                                                                                                           |
| ----------------------------------- | ------------------------------------------------ | ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Introduction                        | Concept and feature overview information         | People interested what REMix does and how we apply it | Summary what REMix is and what can be done with it                                                                                                               |
| Tutorials and example applications  | Installing REMix and learning how to run a model | Users that should be introduced to a functionality    | <ul><li>Python scripts with step-by-step instructions on how to prepare the model</li><li>We decide what the docu reader shall know</li></ul>                    |
| Explanations                        | Understanding concrete REMix concepts            | Advanced users, developers                            | module descriptions: <ul><li>purpose</li><li>required inputs</li><li>constraints and variables</li><li>what can(not) be changed/configured by the user</li></ul> |
| Developing REMix                    | Contributing to REMix development                | Experienced users                                     | <ul><li>Contribution guidelines</li><li>list of contributors</li><li>changelog</li><li>releases</li><ul>                                                         |
| Guides                              | Learning how to do sth.                          | Experienced users                                     | <ul><li>Step-by-step instruction which could be the answer on a user's question</li><li>The user decides what s/he wants to know</li></ul>                     |

### Fixing typos and other smaller mistakes

If you find mistakes are typos in the documentation, you are invited to fix them
if you like.
To do so, you can simply click the "Edit this page" button and follow along.
After the changes are committed to the repository and approved, the
documentation will be rebuilt automatically.


### Generating new content

When wanting to generate new content for the documentation and not just fix
existing one, it might be easier to create a merge request from a GitLab issue.
If you are creating a new branch, make sure to prefix its name with `pages/`, so
that it will automatically be built and appear in the drop-down menu in the top
right corner of the documentation page.
For a general introduction on how to contribute to REMix, head to
{ref}`this section <contributing_active_dev_label>`.

### Adding a new page

To add a new page, add the .md file in the appropriate section (tutorials,
guides, etc.) in the `/docs/` directory and include it in the table of contents.
To include it, simply add the link to the document in one of the tocs of
the `contents.rst` file or any of the files included in the subsections (e.g.
`tutorials/advanced.md`).

Every .md file needs to include a header section.
This includes the title to be displayed as well as the language of the document.
Each page has to start with a level 1 heading (#).

```yaml
---
title: Installing REMix
lang: en-US
---
# Installing REMix
[...]
```

### Building the documentation locally

Instead of having to upload the documentation every time you want to see how
your changes would look like in the final page, it is handy to test them
locally before to save you time and avoid breaking something in the public docu.
To do so, just run

```bash
sphinx-autobuild docs public
```

in the documentation root directory (typically `remix/framework`).

After the command has finished executing, you will get a link
(e.g. http://127.0.0.1:8000) which you can open in your browser to see a local
version of the website.
Any changes you do to your local files in the cloned folder will immediately appear there after saving.

N.B.: The changelog is built automatically in the CI pipeline in GitLab, so will
not appear in your local build.
This is not a bug.

[//]: # (On https://gitlab.dlr.de/remix/framework/-/merge_requests/166, there is a description on how to build it locally anyway.)

### Checking if links are working

Sphinx offers the option to check for broken links.
You can build the docs and check for links using the following command:

```
sphinx-build -b linkcheck docs public
```

### Creating figures in the REMix design

If you happen to need a figure to better explain whatever concept you want to
describe in the documentation, there is a figure template you can build on:
{download}`REMix_figures_template_guidelines.svg </img/REMix_figures_template_guidelines.svg>`

---

The documentation is built with [Sphinx](https://www.sphinx-doc.org/) and
the [PyData Sphinx theme](https://github.com/pydata/pydata-sphinx-theme).

### Referencing within docs

It is possible to define link anchors with

```(link_anchor_label)=```

and reference to these link anchors using

``` {ref}`Text for link <link_anchor_label>` ```.

This has the advantage that independent where the document lies we want to
reference to, we can simply insert the link anchor label.
The links do not break on changes in the directory structure.
Also, the links can be inserted in the middle of each document, so if you want
to specifically reference a subsection, table, figure or equation, you can do
that.

To reference sections in the GAMS code documentation (`remix.framework.model`
code) there are default anchors for every page following this structure:

``` {ref}`Text for link <remix_model_{FEATURE}_label>` ``` where `{FEATURE}` has
to be replaced with `core_transfer`, `core_sets`, etc.
If you want to include new labels for subsections of the documentation, you have
to create a new label inside the docstrings of the GAMS code.

To reference sections in the Python code, a similar syntax can be used, but
instead of `{ref}` at the beginning, you need to use `{func}`, `{class}`,
`{meth}` for Python functions, classes or methods respectively depending on what
you actually want to reference.
A functioning cross reference to a `run_remix` function would be written as
follows, for example:

``` {func}`remix.framework.api.run.run_remix` ```

More information about this is available in the
[Sphinx documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-python-objects).

### Technical documentation guide

The inline code documentation of the model is used to parse the technical
documentation section.
There are various things that have to be taken into account when modifying the
code in order for this to work properly.

#### Markdown in code

Most of the functionality of the documentation comes from having markdown
snippets in the source code.
Everything written after the following pattern will be considered markdown code:

```
* //
```

Inside these markdown snippets one can include {identifier} placeholders which
can be used during the parsing to insert tables or equations.
Example:

```
* // ### converter_capacityParam
* // Capacity parameters describe the expansion boundary conditions.
* // {table_converter_capacityParam}
```

Here, the table_converter_capacityParam placeholder will be replaced with the
table with that identifier.
In the following section we explain how to assign objects to the identifiers.
First, we explain the tables.

#### Table regex patterns

The way in which we find the tables from the code is using
`regular expressions`, also known as `regex`.
The pattern matching for the tables is rather sensitive, it will always look for
the following structure:

```
* // {table_x}
set pc_x
    /
    parameter_1          "Title|Description|constraint_a:value;constraint_b:[fields]"
    parameter_2          "Title|Description|"
    /;
table ...
... ;
```

It is rather specific, but done intentionally that way in order to have as
little redundant information in the code as possible.
For example: we are trying to avoid defining the table information twice, this
would mean having to watch out for multiple instances of the same data.

The contents of the parameter_x lines will be then parsed and assigned to a
table which is going to get the identifier table_x.

If you want to implement a new parameter, as long as you follow this pattern the
documentation should be parsed correctly.

#### Equations

The equations are being parsed directly from the code as well, these get
assigned placeholders equivalent to their definitions.
Notice the previously mentioned placeholder format, example:

```
* // {Eq_converter_unitsBalance}
Eq_converter_unitsBalance(nodesModelSel,years,converter_techs,vintage)
    $(converter_usedTech(nodesModelSel,years,converter_techs,vintage)
        or converter_usedTech(nodesModelSel,years-1,converter_techs,vintage))
    ...
```

To get the equations from the code there is a tool provided directly by GAMS
named "model2tex".
This will turn the equations into LaTeX expressions.
To use it you would first need to run the model for it to produce a
documentation file.

```
remix run --datadir=testing/instances/minimal_lp/data --docfile=docs/scripts/remix --roundts=1
model2tex.sh docs/scripts/remix -m remix
```

We can't use these outputs directly in markdown so they have to be further
processed.
This is done during the parsing step.
The script will take care of assigning the equations to their respective
placeholders in the documentation.

```
python ./docs/scripts/parse_model_markdown.py
```

#### Metadata

The model has in some places some flags that may at first look confusing.
For example:

```
** // SET: accNodes|Source/Accounting nodes|set_accnodes.csv
```

Notice the double asterisk.
This won't be parsed as a markdown file but it is necessary to generate metadata
related to the input files.
If you want to add a new set input file, it is encouraged that you add a similar
pattern for it to be taken into account during validation.
The first field is the name the set would have in a table, the second is a more
human-friendly title and the third is the name of the csv(dat) file you would
find the set on.

There are also some other examples of metadata strings to be highlighted.

```
** // SAME: accLinksData|linksData
```

The SAME flag will assign the first element to the second one at the time of
validating the input data.

```
** // MAP: nodesData|accNodes|map_accnodes.csv
```

The MAP flag will construct a mapping between the first and the second set and
assign it to the given csv(dat) file name.

```
** // UNION: accNodesData|Accounting and data nodes|set_accnodesdata.csv
```

The UNION flag is assigned to the union of one or more input sets.
It won't actually generate a union but it is important that you declare it if
you are using this kind of construction.

#### Notes on regex

The parsing of these elements is rather sensitive to the current model
structure.
So if you are changing some generic structure of the model you better have a
good reason.
Please do a merge request if this is the case.
If it breaks the documentation, you can always try to figure out a solution in
the parsing script.

## Tips on contributing to the source code

...
