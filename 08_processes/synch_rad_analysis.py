import pybdsim
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def analysis(file_name = None) :
    '''
    Plot 2D histogram of energy/angle correlation
    '''

    d = pybdsim.DataPandas.Load(file_name)

    photon_energy = []
    photon_angle  = []

    # loop over events
    print(len(d.get_events()))
    for ievt, ntraj in enumerate(d.get_events()['ntraj']):

        print(ievt)

        for itraj in range(0,ntraj) :

            # get trajectory
            t, td = d.get_trajectory(ievt, itraj)

            # first step
            first_step = t.iloc[0]

            # continue loop if not photon
            if td['parentID'] != 0 :
                photon_angle.append(first_step['py']/first_step['pz'])
                photon_energy.append(first_step['kineticEnergy'])


    photon_energy = numpy.array(photon_energy)
    photon_angle = numpy.array(photon_angle)/1e-3

    xbins = numpy.logspace(numpy.log10(photon_energy.min()),numpy.log10(photon_energy.max()), 50)
    ybins = numpy.linspace(photon_angle.min(), photon_angle.max(), 50)
    counts, xedges, yedges = numpy.histogram2d(photon_energy, photon_angle, bins=[xbins,ybins])  # angle in mrad


    ax = plt.figure(figsize=(7, 7))

    plt.subplot(2,2,1)
    plt.pcolormesh(xedges, yedges, counts.T, cmap='viridis', norm=LogNorm())
    plt.gca().set_xscale('log')
    plt.xlabel('energy/GeV')
    plt.ylabel('angle/mrad')

    plt.subplot(2,2,2)
    plt.hist(photon_angle,50,(photon_angle.min(), photon_angle.max()), orientation='horizontal')

    plt.subplot(2,2,3)
    plt.hist(photon_energy,50,(photon_energy.min(), photon_energy.max()))

    plt.show()
