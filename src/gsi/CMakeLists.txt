cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    GasSurfaceInteraction.cpp
)

install(FILES DrivingForces.h DESTINATION include/mutation++)
install(FILES GasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES GSIReaction.h DESTINATION include/mutation++)
install(FILES SurfaceProductionTerms.h DESTINATION include/mutation++)
install(FILES SurfaceProperties.h DESTINATION include/mutation++)
