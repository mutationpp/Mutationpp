cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    MillikanWhite.cpp
    OmegaVT.cpp
    OmegaET.cpp
    OmegaCV.cpp
    OmegaCElec.cpp
    OmegaCE.cpp
)

install(FILES MillikanWhite.h DESTINATION include/mutation++)
install(FILES TransferModel.h DESTINATION include/mutation++)
