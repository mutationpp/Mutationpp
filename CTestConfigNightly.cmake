##########################################################################
# Site-specific setup
##########################################################################

#execute_process(COMMAND echo "THE GCC OF THIS MACHINE IS EVIRONMENT = $ENV{MY_GCC_VERSION}" )
# Dashboard model (Continuous, Experimental, Nightly)
set(MODEL Continuous)
# source dir(<> dev directory)
set(CTEST_SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src")
# build dir
set(CTEST_BINARY_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/src")

set(MY_TEST_PATH ${CTEST_BINARY_DIRECTORY}/Testing/Temporary)

# commande that containes information about machine and status 
execute_process(COMMAND id -un  COMMAND tr -d "\n" OUTPUT_VARIABLE myName)
set(myMachine ${CMAKE_SYSTEM})
execute_process(COMMAND echo "${myName}@${myMachine}" COMMAND sed "s/ /_/g" OUTPUT_VARIABLE MY_SITE)

#CMAKE_CXX_COMPILER_VERSION
# to get the gcc version 

set(CF_BUILD_TYPE "Release")


execute_process(COMMAND gcc -dumpversion OUTPUT_VARIABLE GCC_VERSION)
execute_process(COMMAND echo "gcc-${MY_GCC_VERSION}" OUTPUT_VARIABLE MY_BUILD)



# Build and site identification
set(CTEST_SITE "${MY_SITE}")
set(CTEST_BUILD_NAME "gcc-${GCC_VERSION}")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")


# commands
set(CTEST_CMAKE_COMMAND "/usr/bin/cmake")
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


 set(START_TIME ${CTEST_ELAPSED_TIME})
ctest_start (${MODEL})
 ctest_update(RETURN_VALUE HAD_UPDATES)
 if(${HAD_UPDATES} GREATER 0)
#  ctest_configure()
#  ctest_build()		
#  ctest_test()
#  ctest_submit()
  execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mailto ${MY_TEST_PATH})#${MY_MAIL_TO} )
 else(${HAD_UPDATES} GREATER 0)
	MESSAGE("There are no new updates since last regression test")	
 endif(${HAD_UPDATES} GREATER 0)

#execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mailto ${MY_TEST_PATH})#${MY_MAIL_TO} )




 if(${HAD_UPDATES} GREATER 0)
 	MESSAGE("There are some new updates since last regression test ==> ${HAD_UPDATES} ") 
 else()
       MESSAGE("There are no new updates since last regression test ==> ${HAD_UPDATES}")  
 endif()



