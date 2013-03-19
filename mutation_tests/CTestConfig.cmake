## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.

set(CTEST_PROJECT_NAME "Mixtures")
set(CTEST_NIGHTLY_START_TIME "23:00:00 GMT")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "localhost")
set(CTEST_DROP_LOCATION "/submit.php?project=mutationtest&subproject=Mixtures")
set(CTEST_DROP_SITE_CDASH TRUE)

set(CTEST_UPDATE_TYPE "svn")
