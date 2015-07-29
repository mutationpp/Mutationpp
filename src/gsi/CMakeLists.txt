cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    CatalysisReaction.cpp
    DiffusionVelocityCalculator.cpp
    GasSurfaceInteraction.cpp
    GSIRateLawGammaConst.cpp
    GSIRateManagerGamma.cpp
    GSIStoichiometryManager.cpp
    SurfaceBalanceSolver.cpp
    SurfaceDescription.cpp
    SurfaceProductionRate.cpp
    SurfacePropertiesGamma.cpp
    WallState.cpp
)

install(FILES DataWallProductionTerms.h DESTINATION include/mutation++)
install(FILES DataGSIRateLaw.h DESTINATION include/mutation++)
install(FILES DataGSIRateManager.h DESTINATION include/mutation++)
install(FILES DataGSIReaction.h DESTINATION include/mutation++)
install(FILES DataGSIRateManager.h DESTINATION include/mutation++)
install(FILES DiffusionDrivingForces.h DESTINATION include/mutation++)
install(FILES DiffusionVelocityCalculator.h DESTINATION include/mutation++)
install(FILES GasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES GSIRateLaw.h DESTINATION include/mutation++)
install(FILES GSIRateManager.h DESTINATION include/mutation++)
install(FILES GSIReaction.h DESTINATION include/mutation++)
install(FILES GSIStoichiometryManager.h DESTINATION include/mutation++)
install(FILES MassBlowingRate.h DESTINATION include/mutation++)
install(FILES SurfaceBalanceSolver.h DESTINATION include/mutation++)
install(FILES SurfaceDescription.h DESTINATION include/mutation++)
install(FILES SurfaceProperties.h DESTINATION include/mutation++)
install(FILES WallProductionTerms.h DESTINATION include/mutation++)
install(FILES WallState.h DESTINATION include/mutation++)
