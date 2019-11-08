from multiplierz.mzReport import ReportReader
from multiplierz.mzTools.mzIdentMLAPI import mzIdentML


class mzIdentMLReader(ReportReader):
    def __init__(self, file_name):
        self.file_name = file_name
        self.datafile = mzIdentML(file_name)
        self.data = self.datafile.peptideSummary()
        self.columns = list(self.data[0].keys())
        
    
    #def __iter__(self):
        #for entry in self.data:
            #vals = entry.values()
            #for i in range(0, len(vals)):
                #if type(vals[i]) == type([]):
                    #vals[i] = "; ".join(vals[i])
                    
            #yield entry.values()
    def __iter__(self):
        for entry in self.data:
            for field in entry:
                if isinstance(entry[field], list):
                    entry[field] = '; '.join(entry[field])
            yield entry
            
    def close(self):
        del self.datafile
        del self.data
        
        

    