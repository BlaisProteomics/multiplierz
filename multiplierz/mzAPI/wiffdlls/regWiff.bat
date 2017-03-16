echo Registering WiffReaderCOM...
@echo off

cd c:\windows\Microsoft.NET\

if exist Framework64 (
    cd Framework64
) ELSE (
    cd Framework
)

cd v4*
if errorlevel 1 (
    cd v3*
    if errorlevel 1 (
    cd v2*    
    )
)

@echo on
RegAsm %1 \tlb \codebase