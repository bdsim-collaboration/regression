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

    counts, xedges, yedges = numpy.histogram2d(photon_energy, photon_angle, bins=50)

    plt.figure(figsize=(7, 6))
    plt.imshow(
        numpy.log10(counts.T),           # transpose: histogram2d returns [x_bin, y_bin], imshow wants [row, col]
        origin='lower',
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        aspect='auto',
        cmap='viridis'
    )
    plt.colorbar(label='Counts')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D Histogram')
    plt.show()

    return photon_energy, photon_angle
