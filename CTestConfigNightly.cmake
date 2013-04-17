##########################################################################
# Site-specific setup
##########################################################################

# Dashboard model (Continuous, Experimental, Nightly)
set(MODEL Continuous)

# source dir(<> dev directory)
set(CTEST_SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src")

# build dir
set(CTEST_BINARY_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/src")

# commande that containes information about machine and status 


execute_process(COMMAND id -un  COMMAND tr -d "\n" OUTPUT_VARIABLE myName)
#execute_process(COMMAND cat /etc/issue.net  COMMAND tr -d "\n" OUTPUT_VARIABLE myMachine)
set(myMachine ${CMAKE_SYSTEM})
execute_process(COMMAND echo "${myName}@${myMachine}" COMMAND sed "s/ /_/g" OUTPUT_VARIABLE MY_SITE)

# to get the gcc version 
set(CF_BUILD_TYPE "Release")
execute_process(COMMAND gcc --version  COMMAND sed -rn "s/gcc[^[:digit:]]*([0-9.]*)/\\1/p"
		COMMAND awk "{print $NF}"  COMMAND  sed -e "s/[^0-9.-]//g" OUTPUT_VARIABLE compiler)
# execute_process(COMMAND echo "gcc-${compiler}-${CF_BUILD_TYPE}-Main" COMMAND sed "s/ /_/g" OUTPUT_VARIABLE MY_BUILD)
execute_process(COMMAND awk "/define BOOST_LIB_VERSION/ {print \$3}" /usr/include/boost/version.hpp
	        COMMAND  sed "s/\"//g" 
		COMMAND	sed "s/_/./g"
		OUTPUT_VARIABLE MY_BOOST)
execute_process(COMMAND echo "gcc-${compiler}-${CF_BUILD_TYPE}-Main" 
		COMMAND sed "s/ /_/g" 
		OUTPUT_VARIABLE MY_COMPILER)

execute_process(COMMAND echo "${MY_COMPILER}  BOOST VERSION ${MY_BOOST}" 
		COMMAND tr -d "\n"
		OUTPUT_VARIABLE MY_BUILD)
# Build and site identification
set(CTEST_SITE "${MY_SITE}")
set(CTEST_BUILD_NAME "${MY_BUILD}")
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


#while (${CTEST_ELAPSED_TIME} LESS 72000)
 set(START_TIME ${CTEST_ELAPSED_TIME})
 ctest_start (${MODEL})
 #ctest_update(RETURN_VALUE HAD_UPDATES)
# if(${HAD_UPDATES} GREATER 0)
  ctest_configure()
  ctest_build()		
  ctest_test()
  ctest_submit()
# endif()
# ctest_sleep( ${START_TIME} 60 ${CTEST_ELAPSED_TIME})
#endwhile()


