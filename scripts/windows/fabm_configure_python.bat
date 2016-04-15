@rem Script to configure the fabm0d executable uing CMake

@set old=%cd%

@echo Build directory:
@if "%build_dir%"=="" ( @set build_dir=%TEMP%\build\fabm ) else ( @echo build_dir is set )
@echo %build_dir%
@if not EXIST "%build_dir%\." ( @mkdir "%build_dir%" )
@chdir "%build_dir%"

@echo Base directories:
@set FABM_BASE=%USERPROFILE%\Documents\FABM\code
@set GOTM_BASE=%USERPROFILE%\Documents\GOTM\code
@echo %FABM_BASE%
@echo %GOTM_BASE%

@echo Default Fortran compiler is ifort
@set compiler=ifort

@echo Install directory:
@echo using python specific install directory

@set host=python
@if not EXIST "%host%\." ( mkdir "%host%" )
@chdir "%host%"
@echo Ready to configure:
cmake "%FABM_BASE%\src\drivers\python" ^
      -DFABM_EMBED_VERSION=on ^
      -DCMAKE_Fortran_COMPILER=%compiler%

@pause

@chdir ..\
@chdir %old%
