import pybdsim
import numpy
import matplotlib.pyplot as plt

def analysis(file_name = None) :
    '''
    Plot 2D histogram of energy/angle correlation
    '''

    d = pybdsim.DataPandas.Load(file_name)

    photon_energy = []
    photon_angle  = []

    # loop over events
    for ievt, ntraj in enumerate(d.get_events()['ntraj']):
        for itraj in range(0,ntraj) :

            # get trajectory
            t = d.get_trajectory(ievt, itraj)

            # first step
            first_step = t.iloc[0]

            # continue loop if not photon
            if first_step['parentID'] != 0 :
                photon_angle.append(first_step['py']/first_step['pz'])
                photon_energy.append(first_step['kineticEnergy'])

    photon_energy = numpy.array(photon_energy)
    photon_angle = numpy.array(photon_angle)/1e-3

    xbins = numpy.logspace(numpy.log10(photon_energy.min()),numpy.log10(photon_energy.max()), 50)
    ybins = numpy.linspace(photon_angle.min(), photon_angle.max(), 50)
    counts, xedges, yedges = numpy.histogram2d(photon_energy, photon_angle, bins=[xbins,ybins])  # angle in mrad

    ax = plt.figure(figsize=(7, 6))
    plt.imshow(
        numpy.log10(counts.T),           # transpose: histogram2d returns [x_bin, y_bin], imshow wants [row, col]
        origin='lower',
        extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
        aspect='auto',
        cmap='viridis'
    )
    plt.colorbar(label='Counts')
    plt.xlabel('energy/GeV')
    plt.ylabel('angle/mrad')
    plt.title('2D Histogram')
    plt.show()

    return photon_energy, photon_angle
