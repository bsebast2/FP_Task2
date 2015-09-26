@echo off
set MATLAB=C:\Program Files\MATLAB\R2015a
"%MATLAB%\bin\win64\gmake" -f Date0501_QuadModel.mk  ISPROTECTINGMODEL=NOTPROTECTING
