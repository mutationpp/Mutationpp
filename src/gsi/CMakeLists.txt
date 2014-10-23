cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    CatalysisRateLaws.cpp
    GasSurfaceInteraction.cpp
    GSIReaction.cpp
)

install(FILES CatalysisRateLaws.h DESTINATION include/mutation++)
install(FILES GasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES GSIReaction.h DESTINATION include/mutation++)
