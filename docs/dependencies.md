<a id="top"></a>

`Mutationpp` depends on `Eigen` and on `Catch2` for testing. Generally speaking the suggested approach is to use these dependencies installing them using the system package manager (or contact your system administrator to make them available to you, for instance as environment module). If they're available in the system in standard locations (or their module is loaded in the current shell), they should be retrieved automatically by CMake (if it does not, it's probably a bug, so open an issue on Github).

If you don't have admin rights on your machine and you cannot get a module, we provide versions in the `thirdparty` directory. You can also upload them using the following procedure:


# Updating Dependencies

The source code of third-party dependencies is located in the `thirdparty` directory. The code for the following dependencies is taken directly from a git repository using the [subrepo](https://github.com/ingydotnet/git-subrepo) utility:

- Eigen ([GitHub mirror of the BitBucket Mercurial repository](https://github.com/RLovelett/eigen))
- Catch ([original GitHub project page](https://github.com/philsquared/Catch))

The installation process of the `subrepo` utility is detailed on the project page.

The status of `subrepo`-managed dependencies can be viewed by issuing the following command from the root of the Mutation++ repository:

```bash
git subrepo status --all
```

## Updating dependency code

If a dependency is managed using `subrepo`, it can be updated by using it. From the root of the Mutation++ repository, issue the following command:

```bash
git subrepo checkout <path_to_dependency> -b <tag_or_branch>
```

For instance, to upgrade Eigen to v3.3.2:

```bash
git subrepo checkout thirdparty/eigen -b 3.3.2
```

Check on the hosting git repo the tag or branch name to use.

If the code is not managed using `subrepo`, delete the current code and replace it with the new version.

## Adding dependency code

If the code is to be managed using `subrepo` (it requires that it is hosted by a git-based service), issue the following command from the root of the Mutation++ repository:

```bash
git subrepo clone <url_to_git_repo> thirdparty/<dependency_directory_name> -b <branch_or_tag>
```

For instance, checking out the Catch source code was done with the following command:
    
```bash
git subrepo clone https://github.com/philsquared/Catch.git thirdparty/catch -b v1.9.3
```

## When using `subrepo`

The `subrepo` utility adds remotes to your working repository. They should be removed when you are done by issuing the following command from the root of the Mutation++ repository:

```bash
git subrepo clean --all
```

If some tags persist, remove all of them with the following command (the Mutation++ tags will be restored next time you fetch from the central repository):

```bash
git tag | xargs git tag -d
```

## Going further

- [More on Git subtrees](https://medium.com/@porteneuve/mastering-git-subtrees-943d29a798ec)
