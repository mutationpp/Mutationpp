cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    MillikanWhite.cpp
    TransferModelVT.cpp
    TransferModelET.cpp
    TransferModelEV.cpp
    TransferModelCV.cpp
)

install(FILES MillikanWhite.h DESTINATION include/mutation++)
install(FILES TransferModel.h DESTINATION include/mutation++)
