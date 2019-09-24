@echo off
REM automatically generated
ECHO generated on host: wpia-didelt319.dide.ic.ac.uk
ECHO generated on date: 2019-09-23
ECHO didehpc version: 0.2.9
ECHO context version: 0.1.4
ECHO running on: %COMPUTERNAME%
set CONTEXT_WORKDRIVE=V:
set CONTEXT_WORKDIR=os210
set CONTEXT_ROOT=V:\os210\contexts
set CONTEXT_ID=4b99c4c8cab6437cf031dd0aece61ea9
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
ECHO mapping V: -^> \\fi--san03.dide.ic.ac.uk\HOMES
net use V: \\fi--san03.dide.ic.ac.uk\HOMES /y
ECHO This is a parallel job: will use %CPP_NUMCPUS%
set CONTEXT_CORES=%CCP_NUMCPUS%
set REDIS_HOST=12.0.0.1
set REDIS_URL=redis://12.0.0.1:6379
%CONTEXT_WORKDRIVE%
cd \%CONTEXT_WORKDIR%
ECHO working directory: %CD%
ECHO this is a single task
set CONTEXT_TASK_ID=8ad240097f958d22282889b24bbefc0d
set CONTEXT_LOGFILE=V:\os210\contexts\logs\%CONTEXT_TASK_ID%
ECHO logfile: %CONTEXT_LOGFILE%
@REM The quoting here is necessary for paths with spaces.
ECHO on
Rscript "V:\os210\contexts\bin\task_run" "%CONTEXT_ROOT%" %CONTEXT_TASK_ID% > "%CONTEXT_LOGFILE%" 2>&1
@ECHO off
%SystemDrive%
set ErrorCode=%ERRORLEVEL%
ECHO Removing mapping Q: -^> \\fi--san03.dide.ic.ac.uk\homes\os210
net use Q: /delete /y
ECHO Removing mapping V: -^> \\fi--san03.dide.ic.ac.uk\HOMES
net use V: /delete /y
set ERRORLEVEL=%ErrorCode%
if %ERRORLEVEL% neq 0 (
  ECHO Error running task
  EXIT /b %ERRORLEVEL%
)
@ECHO Quitting
