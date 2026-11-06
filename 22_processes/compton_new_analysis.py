import pybdsim
import numpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

class AnalysisClass(pybdsim.RootEventAnalysis.RootEventAnalyser) :
    def __init__(self) :
        super().__init__()

    def init(self):
        self.eleX = []
        self.eleY = []
        self.phoX = []
        self.phoY = []
        self.histograms    = {}

    def process(self) :
        sampler = self.event.GetSampler(b'sep.')
        for i in range(sampler.n):
            if sampler.partID[i] == 11:
                self.eleX.append(sampler.x[i]*1e3) # Convert to mm
                self.eleY.append(sampler.y[i]*1e3)
            elif sampler.partID[i] == 22:
                self.phoX.append(sampler.x[i]*1e3)
                self.phoY.append(sampler.y[i]*1e3)

    def terminate(self):
        self.eleX = numpy.array(self.eleX)
        self.eleY = numpy.array(self.eleY)
        self.phoX = numpy.array(self.phoX)
        self.phoY = numpy.array(self.phoY)

        ele_xbins = numpy.linspace(self.eleX.min(),self.eleX.max(), 100)
        ele_ybins = numpy.linspace(self.eleY.min(),self.eleY.max(), 100)
        pho_xbins = numpy.linspace(self.phoX.min(),self.phoX.max(), 100)
        pho_ybins = numpy.linspace(self.phoY.min(),self.phoY.max(), 100)

        ele_xcounts, ele_xedges = numpy.histogram(self.eleX, bins=ele_xbins)
        ele_ycounts, ele_yedges = numpy.histogram(self.eleY, bins=ele_ybins)
        pho_xcounts, pho_xedges = numpy.histogram(self.phoX, bins=pho_xbins)
        pho_ycounts, pho_yedges = numpy.histogram(self.phoY, bins=pho_ybins)

        self.histograms["electron_xprofile"] = {"type": "histo1d",
                                                "contents": ele_xcounts,
                                                "xedges": ele_xedges,
                                                "xlabel": "X [mm]",
                                                "ylabel": "Counts"}

        self.histograms["electron_yprofile"] = {"type": "histo1d",
                                                "contents": ele_ycounts,
                                                "xedges": ele_yedges,
                                                "xlabel": "Y [mm]",
                                                "ylabel": "Counts"}

        self.histograms["photon_xprofile"] = {"type": "histo1d",
                                              "contents": pho_xcounts,
                                              "xedges": pho_xedges,
                                              "xlabel": "X [mm]",
                                              "ylabel": "Counts"}

        self.histograms["photon_yprofile"] = {"type": "histo1d",
                                              "contents": pho_ycounts,
                                              "xedges": pho_yedges,
                                              "xlabel": "Y [mm]",
                                              "ylabel": "Counts"}

        self.add_persistent_data("histograms",self.histograms)


def analysis(file_name  = None) :
    a = pybdsim.RootEventAnalysis.RootEventAnalysis(file_name)
    ac = AnalysisClass()
    a.analysis(ac)
    return ac
