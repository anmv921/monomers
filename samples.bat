@echo off

rem Numero de amostras
SET /A n = 1

rem Ler as sementes
set /A j=0
for /f %%a in (seeds.txt) do (
	set /A j+=1
	call set array[%%j%%]=%%a
	call set m=%%j%%
)

set list=121 144 169 196 225 256 289 324 361
(for %%N in (%list%) do ( 
	echo %%N 
	for /l %%x in (1, 1, %n%) do (
		echo %%x

		call set folder=data\N\%%N\samples\%%x\
		if not exist %%folder%% call mkdir %%folder%%

		call copy params1.in %%folder%%params.in
		call copy main.exe %%folder%%

		call cd %%folder%%
		call echo seed	%%array[%%x]%%>>params.in
		call echo NCell	%%N>>params.in
		main.exe
		cd ..\..\..\..\..
	)  
))

