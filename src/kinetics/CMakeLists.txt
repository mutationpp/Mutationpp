cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    JacobianManager.cpp
    Kinetics.cpp
    MillikanWhite.cpp
    RateLaws.cpp
    RateManager.cpp
    Reaction.cpp
    StoichiometryManager.cpp
)

install(FILES JacobianManager.h DESTINATION include/mutation++)
install(FILES Kinetics.h DESTINATION include/mutation++)
install(FILES MillikanWhite.h DESTINATION include/mutation++)
install(FILES RateLaws.h DESTINATION include/mutation++)
install(FILES RateManager.h DESTINATION include/mutation++)
install(FILES Reaction.h DESTINATION include/mutation++)
install(FILES StoichiometryManager.h DESTINATION include/mutation++)
install(FILES ThirdBodyManager.h DESTINATION include/mutation++)

