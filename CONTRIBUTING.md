Contribution Guide
==================

First of all, thank you for your enthusiasm at contributing to Mutation++. Please 

## How can I help?

You can contribute to Mutation++ in a variety of ways:

* Reporting issues
* Suggesting new features
* Proposing fixes and code improvements

## Workflow in a nutshell

Choice has been made to follow a workflow close to the [GitHub flow](https://guides.github.com/introduction/flow/index.html). This simple workflow displays the following features:

* Heavy use of branching: any new feature of fix is implemented in a dedicated branch, which is then merged into the `master` branch.
* Limited direct commits to the central repository: all contributors work on their [own fork](https://sync.vki.ac.be/help/workflow/forking_workflow.md) of the project and submit their changes through merge requests.

The central repository has a minimal number of active branches, the reference being the `master` branch. This workflow might evolve into something more similar to the [GitLab flow](https://about.gitlab.com/2014/09/29/gitlab-flow/) in the future, depending on whether or not the need for release tracking emerges. The [feature branch workflow](https://sync.vki.ac.be/help/workflow/workflow.md) summarizes the process.

### Getting started

If you want to contribute to Mutation++, fork the central repository and branch the `master` branch to a feature branch. To keep your fork in sync with the upstream repository, define an `upstream` remote and regularly update your repository (see [these instructions](https://help.github.com/articles/syncing-a-fork/)).

### Submitting modifications to the central repository

When you deem it is time to have your changes merged back into the central repository, submit a merge request. Your code will be reviewed and you will get feedback on your work. Code which fails to pass tests will not be merged.

## Coding style

### C++ coding style

Your committed work must adhere to the rules defined in the C++ Coding Style guide.
