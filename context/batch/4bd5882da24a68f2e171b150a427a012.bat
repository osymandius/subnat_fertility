@echo off
REM automatically generated
ECHO generated on host: WPIA-DIDE215
ECHO generated on date: 2020-03-09
ECHO didehpc version: 0.2.9
ECHO context version: 0.1.4
ECHO running on: %COMPUTERNAME%
set CONTEXT_WORKDRIVE=Q:
set CONTEXT_WORKDIR=
set CONTEXT_ROOT=Q:\context
set CONTEXT_ID=d21659cf13544dca433d95551fd06e13
set CONTEXT_PROPAGATE_ERROR=TRUE
set CONTEXT_BOOTSTRAP=TRUE
call setr64_3_6_1.bat
REM If Java is wanted, then call setJava64.
REM If called with blank, it adds default JRE.
IF 'FALSE'=='TRUE' (
  call setJava64.bat 
)
ECHO mapping Q: -^> \\fi--san03.dide.ic.ac.uk\homes\os210
net use Q: \\fi--san03.dide.ic.ac.uk\homes\os210 /y
ECHO This is a parallel job: will use %CPP_NUMCPUS%
set CONTEXT_CORES=%CCP_NUMCPUS%
set REDIS_HOST=12.0.0.1
set REDIS_URL=redis://12.0.0.1:6379
%CONTEXT_WORKDRIVE%
cd \%CONTEXT_WORKDIR%
ECHO working directory: %CD%
ECHO this is a single task
set CONTEXT_TASK_ID=4bd5882da24a68f2e171b150a427a012
set CONTEXT_LOGFILE=Q:\context\logs\%CONTEXT_TASK_ID%
ECHO logfile: %CONTEXT_LOGFILE%
@REM The quoting here is necessary for paths with spaces.
ECHO on
Rscript "Q:\context\bin\task_run" "%CONTEXT_ROOT%" %CONTEXT_TASK_ID% > "%CONTEXT_LOGFILE%" 2>&1
@ECHO off
%SystemDrive%
set ErrorCode=%ERRORLEVEL%
ECHO Removing mapping Q: -^> \\fi--san03.dide.ic.ac.uk\homes\os210
net use Q: /delete /y
set ERRORLEVEL=%ErrorCode%
if %ERRORLEVEL% neq 0 (
  ECHO Error running task
  EXIT /b %ERRORLEVEL%
)
@ECHO Quitting
