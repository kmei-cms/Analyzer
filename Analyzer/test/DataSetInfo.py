import ROOT

class DataSetInfo:
    def __init__(self, basedir, fileName, label, sys=None, processName=None, process=None, rate=None, lumiSys=None):
        self.basedir = basedir
        self.fileName = fileName
        self.label = label
        try: 
            self.file = ROOT.TFile.Open(basedir+fileName)
        except:
            pass
        self.sys = sys
        self.processName = processName
        self.process = process
        self.rate = rate
        self.lumiSys = lumiSys

    def getHisto(self, name):
        return self.file.Get(name)

    def getFile(self):
        return self.file

    def __del__(self):
        self.file.Close()
