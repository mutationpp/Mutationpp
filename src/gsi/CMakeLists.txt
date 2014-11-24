cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    CatalysisRateLaws.cpp
    CatalysisRateManager.cpp
    CatalyticFRCWallProperties.cpp
    GasSurfaceInteraction.cpp
    GSIReaction.cpp
)

install(FILES CatalysisRateLaws.h DESTINATION include/mutation++)
install(FILES CatalysisRateManager.h DESTINATION include/mutation++)
install(FILES CatalyticFRCWallProperties.h DESTINATION include/mutation++)
install(FILES GasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES GSIReaction.h DESTINATION include/mutation++)
