add_executable(test_qr test_qr.cpp)
target_link_libraries(test_qr ${Boost_LIBRARIES})
add_test(QR test_qr)