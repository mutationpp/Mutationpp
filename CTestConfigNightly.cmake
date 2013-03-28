##########################################################################
# Site-specific setup
##########################################################################

# Dashboard model (Continuous, Experimental, Nightly)
set(MODEL Continuous)

# source dir(<> dev directory)
set(CTEST_SOURCE_DIRECTORY "/home/didi/mutation++/branches/dinesh/src")

# build dir
set(CTEST_BINARY_DIRECTORY "/home/didi/mutation++/branches/dinesh/build/src/tests")

# Build and site identification
set(CTEST_SITE "Dinesh-VM-Ubuntu12.04")
set(CTEST_BUILD_NAME "gcc4.6.3-Release-Main")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CF_BUILD_TYPE "Release")

# commands
set(CTEST_CMAKE_COMMAND "/usr/local/bin/cmake")
set(CTEST_UPDATE_COMMAND "/usr/bin/svn")


##########################################################################
# Configuration
##########################################################################

# Initial cache entries
set(CF_CONFIG_OPTIONS " -DCMAKE_BUILD_TYPE:STRING=${CF_BUILD_TYPE}")
#set(CF_CONFIG_OPTIONS "${CF_CONFIG_OPTIONS} -DSITE:STRING=${CTEST_SITE} -DBUILDNAME:STRING=${CTEST_BUILD_NAME}")

# Assemble cmake command line
set(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} ${CF_CONFIG_OPTIONS}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

##########################################################################
# Run the dashboard
##########################################################################

#ctest_start (${MODEL})			# demarrage
#ctest_update(RETURN_VALUE HAD_UPDATES)	# svn update
#ctest_configure()			# cmake (...)
#ctest_build()				# make
#ctest_test()				# ctest = make test
#ctest_submit()				# envoyer au serveur


#while (${CTEST_ELAPSED_TIME} LESS 72000)
 # set(START_TIME ${CTEST_ELAPSED_TIME})
 ctest_start (${MODEL})
 ctest_update(RETURN_VALUE HAD_UPDATES)
# if(${HAD_UPDATES} GREATER 0)
  ctest_configure()
  ctest_build()
  ctest_test()
  ctest_submit()
# endif()
# ctest_sleep( ${START_TIME} 60 ${CTEST_ELAPSED_TIME})
#endwhile()


