# Build and install executables
set(mutation++_exes checkmix mppequil bprime mppshock)
foreach (exe ${mutation++_exes})
    add_executable(${exe} ${exe}.cpp)
    target_link_libraries(${exe}
        PRIVATE
            mutation++
    )
    install(
        TARGETS ${exe}
        EXPORT mutation++Targets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
)
endforeach(exe ${mutation++_exes})
