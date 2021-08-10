from multiplierz.mzReport import ReportReader
from multiplierz.mzTools.mzIdentMLAPI import mzIdentML


class mzIdentMLReader(ReportReader):
    def __init__(self, file_name):
        self.file_name = file_name
        self.datafile = mzIdentML(file_name)
        self.data = self.datafile.peptideSummary()
        try:
            self.columns = list(self.data[0].keys())
        except IndexError:
            self.columns = None # There's no results so this shouldn't case problems.
    
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


if __name__ == '__main__':
    # foo = mzIdentMLReader('C:\\Users\\WilliamMaxAlexander\\Projects\\K562_vs_HUVEC\\K562\\SG_2\\searchgui_out\\20100614_Velos1_TaGe_SA_K562_5.raw.msgf.mzid.gz')
    foo = mzIdentMLReader('/mnt/c/users/WilliamMaxAlexander/projects/K562_vs_HUVEC/K562/SG_2/searchgui_out/20100614_Velos1_TaGe_SA_K562_5.raw.msgf.mzid.gz')
    raise Exception
    print(list(foo))
