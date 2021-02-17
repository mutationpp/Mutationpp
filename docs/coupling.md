<a id="top"></a>

# How to use `mutation++` from your code

`mutation++` uses `CMake` as build toolchain. In particular, starting from
version `0.3.0`, `CMake` exports a `mutation++Config.cmake` file together with
the specific configuration files for the build type (e.g. `Release`, `Debug`)
etc... that allows you to easily include the library into your own code if you
also use `CMake` (ref: [`find_package(mutation++ CONFIG)`][1]). Our examples,
that you can find in the `examples` directory, can be built independently from
`mutation++`. You can take inspiration from the `CMakeLists.txt` that you find
in each example directory in order to require and use `mutation++` in your
code. 

https://github.com/mutationpp/Mutationpp/blob/1d34aec5df89ae7c8d2b76a816e95b02efb4e8bb/examples/c%2B%2B/CMakeLists.txt.shared#L1-L20

So, assuming that you (or your sysadmin) installed `mutation++` in a directory
that is also in the install prefix where `CMake` looks for `*.cmake` files
(ref: [CMAKE_PREFIX_PATH][2]), you basically just need to
`find_package(mutation++)` and then link your target against `mutation++`
target:

```
target_link_libraries(<my-target>
    [PRIVATE|PUBLIC|INTERFACE]
        mutation++::mutation++
)
```


[1]: https://cmake.org/cmake/help/latest/command/find_package.html#full-signature-and-config-mode
[2]: https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html
