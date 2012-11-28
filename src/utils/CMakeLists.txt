cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    Units.cpp
    utilities.cpp
    XMLite.cpp
)

install(FILES AutoRegistration.h DESTINATION include/mutation++)
install(FILES IteratorWrapper.h DESTINATION include/mutation++)
install(FILES LookupTable.h DESTINATION include/mutation++)
install(FILES ReferenceServer.h DESTINATION include/mutation++)
install(FILES Units.h DESTINATION include/mutation++)
install(FILES XMLite.h DESTINATION include/mutation++)

