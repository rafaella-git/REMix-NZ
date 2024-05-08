# REMix tutorials

This repository contains tutorials on building basic models in the energy system modeling framework REMix. This is done with Python (although the framework itself is written in GAMS).

To get an overview and quick introduction of the tutorial content, head over to https://dlr-ve.gitlab.io/esy/remix/framework/dev/getting-started/tutorials/ for the complete guides.

## Additional notes on setup of tutorials

To be able to do the tutorials, you need to have a working REMix installation ready to run on your machine. Make sure you follow the installation instructions on https://dlr-ve.gitlab.io/esy/remix/framework/dev/getting-started/install-remix.html first.

For all tutorials to run correctly, the environment you use needs to fulfill the requirements listed in `tutorial_env.yaml`. If for any reason you do not want to extend the environment created when following the installation instruction (`remix-env`), you can create a dedicated environment from that file e.g. using the mamba package manager as follows (that will be called `remixtutorials`):

```
mamba env create -f tutorial_env.yaml
```

Usually—i.e. (1) if you cloned the REMix repository and did not change file paths afterwards; and (2) if you set your system environment variables for GAMS according to the instructions on https://dlr-ve.gitlab.io/esy/remix/framework/dev/getting-started/install-remix.html#adding-gams-as-environment-variable, the REMix tutorials should run as they are shipped on your machine.

Depending on the IDE of your choice, you might have to change the settings so that scripts are executed in the directory they are in. Otherwise, the relative paths set will not work for you and you will have to adjust them.
