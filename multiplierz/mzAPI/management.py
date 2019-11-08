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
    print("\n-")
    print("Press enter to continue.")
    raw_input()
    return

import sys
#foo = sys.path + [os.path.dirname(sys.executable)]

#print foo
#print zip(*os.walk(foo[1]))[0]
#print os.path.dirname(__file__)
#print sys.executable

msfilereader_installer = os.path.join(os.path.dirname(__file__), 'MSFileReader_x86_x64_v3.0SP3.exe')

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
    print("Checking for access to MS data formats...")
    
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
                        "BlaisT2D.dll" : "T2D",
                        "RawReader.dll" : "RAW"}
    interfaceGUIDs = {"T2D" : "{7e3450b1-75e7-49b2-9be7-64cbb2458c56}",
                      "D" : 'Agilent.MassSpectrometry.DataAnalysis.MassSpecDataReader',
                      "WIFF" : "{9eabbbb3-5a2a-4f73-aa60-f87b736d3476}",
                      "RAW" : "{10729396-43ee-49e5-aa07-85f02292ac70}"}
    
    
    if sys.maxsize > 2**32: bitness = "64"
    else: bitness = ""
    
    try:
        netDir = sorted(glob.glob('c:/Windows/Microsoft.NET/Framework%s/v[34]*/RegAsm.exe' % bitness))[-1]
    except IndexError:
        message =  (".NET version 3 or higher not found; the free .NET redistributable can be found at: "
                    "http://www.microsoft.com/en-us/download/details.aspx?id=17718")
        print(message)
        pause()
        return 1   
    
    initialChecks = {}
    for filetype, guid in list(interfaceGUIDs.items()):
        try:
            Dispatch(guid)
            initialChecks[filetype] = True
        except:
            initialChecks[filetype] = False

    
    dllsFound = {}
    for directory in sys.path + [os.path.dirname(sys.executable)]:
        if len(os.path.abspath(directory)) <= 3: # Searching from base directory is unwise.
            print(("Not searching from %s" % directory))
            continue 
        for path, names, files in os.walk(directory):
            #print path
            for filename in files:
                if (filename in list(interfaceModules.keys())) and (not filename in list(dllsFound.keys())):
                    dllsFound[filename] = os.path.join(path, filename)
                    if len(dllsFound) == len(interfaceModules): break
    
    registerResults = defaultdict(int)
    for filename in list(interfaceModules.keys()):
        print(('\n\n' + filename + '\n'))
        try:
            dllpath = dllsFound[filename]
            ret = subprocess.call([netDir, dllpath, "/tlb", "/codebase"])
            registerResults[interfaceModules[filename]] |= ret
        except KeyError:
            registerResults[interfaceModules[filename]] = "No .DLL"
    
        
    afterChecks = {}
    for filetype, guid in list(interfaceGUIDs.items()):
        try:
            Dispatch(guid)
            afterChecks[filetype] = True
        except:
            afterChecks[filetype] = False
            
    
    print("Registration operations completed.")
    print("Results: \n")
    print("File Type\tRegistered Before\tRegAsm Return Code\tRegistered Now")
    for filetype in list(interfaceGUIDs.keys()):
        print(("{0}\t\t{1}\t\t{2}\t\t\t{3}".format(filetype,
                                                  initialChecks[filetype],
                                                  registerResults[filetype],
                                                  afterChecks[filetype])))
    print("\n")
    #if not afterChecks["RAW"]:
        #print "MSFileReader (required for RAW file access) has not been installed.  Run the Thermo MSFileReader installation package now? [Y/n]"
        #if 'n' not in raw_input().lower():
            #if not os.path.exists(msfilereader_installer):
                #print "MSFileReader installer not found!  Please re-install multiplierz or download MSFileReader from the Thermo Scientific website."
                #print msfilereader_installer
                #return 1
            #print "Please wait..."
            #retval = subprocess.call(msfilereader_installer)
            #if not retval:
                #print "Done."
            #else:
                #print "An error occurred. (Return value %s)" % retval
    
        
    
    return 0

if __name__ == '__main__':
    import ctypes, sys
    import win32com.shell.shell as shell
    from win32com.shell import shellcon
    import win32event, win32process, win32con # con not com?!
    from time import sleep

    # Much Windows interface voodoo;
    # Calls registerInterfaces(), with the complication that if 
    # this is not being run with administrator priviledges, it runs
    # a new Python process calling this module that requests
    # admin priviledges on startup.

    if ctypes.windll.shell32.IsUserAnAdmin():
        registerInterfaces()
        pause()
    else:
        print((sys.executable, sys.argv))
        #reRun = ' '.join([sys.executable] + sys.argv)
        #foo = ctypes.windll.shell32.ShellExecuteW(None, "runas", sys.executable, "", None, 1)
        #procinfo = shell.ShellExecuteEx(lpVerb = 'runas', lpFile = sys.executable,
                                        #lpParameters = ' '.join(sys.argv), nShow = 5 )
        procinfo = shell.ShellExecuteEx(nShow = win32con.SW_SHOWNORMAL,
                                        fMask = shellcon.SEE_MASK_NOCLOSEPROCESS,
                                        lpVerb = 'runas',
                                        lpFile = sys.executable,
                                        lpParameters = ' '.join(sys.argv))
        sleep(0.5)
        obj = win32event.WaitForSingleObject(procinfo['hProcess'], win32event.INFINITE)
        try:
            retval = win32process.GetExitCodeProcess(procinfo['hProcess'])
        except:
            retval = 'No return value.'
        
        print(retval)
        
    