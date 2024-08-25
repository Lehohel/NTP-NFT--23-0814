set list=stRK4 mtRK4, stABM, mtABM, cuRK4, cuABM

(for %%a in (%list%) do (
    cd ./%%a
    make clean build
    .\bin\%%a.exe -output=output_%%a
    cd ..
))

