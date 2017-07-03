# Copyright 2008 Dana-Farber Cancer Institute
# multiplierz is distributed under the terms of the GNU Lesser General Public License
#
# This file is part of multiplierz.
#
# multiplierz is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiplierz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiplierz.  If not, see <http://www.gnu.org/licenses/>.


# Utility commands to make getting multiplierz to work a bit easier.

import os

def pause():
    print "\n-"
    print "Press enter to continue."
    raw_input()
    return

import sys
#foo = sys.path + [os.path.dirname(sys.executable)]

#print foo
#print zip(*os.walk(foo[1]))[0]
#print os.path.dirname(__file__)
#print sys.executable

def testInterfaces():
    from win32com.client import Dispatch
    guids = ['MSFileReader.XRawfile',
             '{9eabbbb3-5a2a-4f73-aa60-f87b736d3476}',
             '{7e3450b1-75e7-49b2-9be7-64cbb2458c56}',
             'Agilent.MassSpectrometry.DataAnalysis.MassSpecDataReader']
    worked = []
    for guid in guids:
        try:
            Dispatch(guid)
            worked.append(True)
        except:
            worked.append(False)
    return worked

def registerInterfaces():
    print "Checking for access to MS data formats..."
    
    import ctypes
    assert ctypes.windll.shell32.IsUserAnAdmin(), "registerInterfaces() must be run with administrator priviledges!"
    
    import os
    import sys
    import subprocess
    import glob
    from collections import defaultdict
    from win32com.client import Dispatch
    
    interfaceModules = {"BlaisWiff.dll" : "WIFF",
                        "MassSpecDataReader.dll" : "D", # Agilent
                        "BaseCommon.dll" : "D",
                        "BaseDataAccess.dll" : "D",
                        "BlaisT2D.dll" : "T2D"}
    interfaceGUIDs = {"T2D" : "{7e3450b1-75e7-49b2-9be7-64cbb2458c56}",
                      "D" : 'Agilent.MassSpectrometry.DataAnalysis.MassSpecDataReader',
                      "WIFF" : "{9eabbbb3-5a2a-4f73-aa60-f87b736d3476}",
                      "RAW" : "MSFileReader.XRawfile"}
    
    
    if sys.maxsize > 2**32: bitness = "64"
    else: bitness = ""
    
    try:
        netDir = sorted(glob.glob('c:/Windows/Microsoft.NET/Framework%s/v[34]*/RegAsm.exe' % bitness))[-1]
    except IndexError:
        message =  (".NET version 3 or higher not found; the free .NET redistributable can be found at: "
                    "http://www.microsoft.com/en-us/download/details.aspx?id=17718")
        print message
        pause()
        return 1   
    
    initialChecks = {}
    for filetype, guid in interfaceGUIDs.items():
        try:
            Dispatch(guid)
            initialChecks[filetype] = True
        except:
            initialChecks[filetype] = False

    
    dllsFound = {}
    for directory in sys.path + [os.path.dirname(sys.executable)]:
        if len(os.path.abspath(directory)) <= 3: # Searching from base directory is unwise.
            print "Not searching from %s" % directory
            continue 
        for path, names, files in os.walk(directory):
            #print path
            for filename in files:
                if (filename in interfaceModules.keys()) and (not filename in dllsFound.keys()):
                    dllsFound[filename] = os.path.join(path, filename)
                    if len(dllsFound) == len(interfaceModules): break
    
    registerResults = defaultdict(int)
    for filename in interfaceModules.keys():
        print '\n\n' + filename + '\n'
        try:
            dllpath = dllsFound[filename]
            ret = subprocess.call([netDir, dllpath, "/tlb", "/codebase"])
            registerResults[interfaceModules[filename]] |= ret
        except KeyError:
            registerResults[interfaceModules[filename]] = "No .DLL"

    registerResults["RAW"] = "N/A"
        
    afterChecks = {}
    for filetype, guid in interfaceGUIDs.items():
        try:
            Dispatch(guid)
            afterChecks[filetype] = True
        except:
            afterChecks[filetype] = False
            
    
    print "Registration operations completed."
    print "Results: \n"
    print "File Type\tRegistered Before\tRegAsm Return Code\Registered Now"
    for filetype in interfaceGUIDs.keys():
        print "{0}\t\t{1}\t\t{2}\t\t\t{3}".format(filetype,
                                                  initialChecks[filetype],
                                                  registerResults[filetype],
                                                  afterChecks[filetype])
    print "\n"
    if not afterChecks["RAW"]:
        print "\n To register RAW file access, download and install the free Thermo MS File Reader."
        return 1
    
    return 0

if __name__ == '__main__':
    registerInterfaces()