set REGASMPROG=%SystemRoot%\Microsoft.NET\Framework64\v4.0.30319\regasm.exe

%REGASMPROG% BaseCommon.dll /tlb /u
%REGASMPROG% BaseDataAccess.dll /tlb /u
%REGASMPROG% MassSpecDataReader.dll /tlb /u

%REGASMPROG% BaseCommon.dll /tlb 
%REGASMPROG% BaseDataAccess.dll /tlb 
%REGASMPROG% MassSpecDataReader.dll /tlb 

