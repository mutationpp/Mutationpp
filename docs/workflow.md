<a id="top"></a>
# Workflow Guide

The Mutation++ development workflow follows the [GitHub flow](https://guides.github.com/introduction/flow/index.html) closely. This simple workflow has the following features:

* Heavy use of branching: any new feature or fix is implemented in a dedicated branch, which is then merged into the `master` branch.
* Limited direct commits to the central repository: all contributors work on their [own fork](https://sync.vki.ac.be/help/workflow/forking_workflow.md) of the project and submit their changes through merge requests.

The central repository has a minimal number of active branches, the reference being the `master` branch. This workflow might evolve into something more similar to the [GitLab flow](https://about.gitlab.com/2014/09/29/gitlab-flow/) in the future, depending on whether or not the need for release tracking emerges. The [feature branch workflow](https://sync.vki.ac.be/help/workflow/workflow.md) summarizes the process.

## Getting Started

If you want to contribute to Mutation++, fork the central repository and [branch the `master` branch to a feature branch](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-new-branch). Branches created to address issues are commonly named after the issue ID and title branch:

    git checkout master
    git checkout -b {issue number}_{branch name}

Example: `29_error_handling_framework` is created to solve issue #29 "Add proper error handling through C++ exceptions").

To keep your fork in sync with the upstream repository, define an `upstream` remote and regularly update your repository (see [these instructions](https://help.github.com/articles/syncing-a-fork/)).

## Submitting Modifications to the Central Repository

When you deem it is time to have your changes merged back into the central repository, submit a [pull request](https://help.github.com/articles/about-pull-requests/). Your code will be reviewed and you will get feedback on your work. Code which fails to pass tests will not be merged.
