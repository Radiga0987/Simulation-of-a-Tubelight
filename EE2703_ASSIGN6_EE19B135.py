#!/usr/bin/env python3
"""The above line is for portability across different systems in case the
language interpreter is installed in different locations."""

# Importing required libraries
import numpy as np
import sys
import matplotlib.pyplot as plt


def Tubelight_sim(n, M, nk, u0, p, Msig):  # Function for simulating the tubelight
    # Initializing the required lists to list of zeros
    xx = np.zeros(n*M)
    u = np.zeros(n*M)
    dx = np.zeros(n*M)

    """The following lists are kept empty initially as
    we do not know their actual lengths.As and when required,
    we append the values to these lists"""
    I = []
    X = []
    V = []

    # Looping over the turns
    for k in range(nk):

        """For the first iteration we check where xx>0.
        But obviously this is not required as locations are
        zero for first iteration.(So this can be removed,but
        I have kept it just to see the logic).For the remaining
        iterations,we dont do this because it has already been
        done at the ending of previous iteration."""
        if k == 0:
            ii = np.where(xx > 0)[0]

        # Moving the electrons that exist by calculating displacement and adding to current position
        dx[ii] = u[ii]+0.5
        xx[ii] += dx[ii]

        # Increasing the velocity of the electrons that exist.
        u[ii] += 1

        # Determining which particles have hit the anode
        anode_hit = np.where(xx >= n)
        # The positions displacements and velocities of these are set to 0
        xx[anode_hit] = 0
        u[anode_hit] = 0
        dx[anode_hit] = 0

        # Determining electrons with velocity more than threshold velocity
        kk = np.where(u >= u0)[0]
        # Of these, which electrons are ionized (Probabilistic process)
        ll = np.where(np.random.rand(len(kk)) <= p)[0]
        kl = kk[ll]

        # Reset the velocities of these electrons(ionized) to zero
        u[kl] = 0
        # The collision could have occurred at any point between the previous xi and the current xi
        xx[kl] -= dx[kl]*np.random.rand()

        """Actually the positions should not be distributed uniformly.
        Electrons are accelerating, and uniformly distributed in
        time does not mean uniformly distributed in space.Better 
        code written below.
        """
        """
        dt = np.random.rand(len(kl))
        xx[kl]=xx[kl]-dx[kl]+((u[kl]-1)*dt+0.5*dt*dt)
        u[kl]=0
        """

        # Excited atoms at this location resulted in emission from that point.
        I.extend(xx[kl].tolist())

        # Inject M new electrons
        m = int(np.random.rand()*Msig+M)

        # Add them to unused slots.
        Slots_add = np.where(xx == 0)[0]

        """We dont need to check if m>len(Slots_add) because the 
        below code will automatically add a maximum of only len(Slots_add)
        number of electrons even if m>len(Slots_add)"""
        xx[Slots_add[0:m]] = 1
        u[Slots_add[0:m]] = 0

        # Finding all the existing electrons.We add their positions and velocities to the X and V vectors.
        existing_electrons = np.where(xx > 0)[0]
        X.extend(xx[existing_electrons].tolist())
        V.extend(u[existing_electrons].tolist())
        # We set ii directly equal to existing electrons here so that we dont have to do it twice every loop
        # This ii will serve as the ii for next iteration
        ii = existing_electrons

    return X, V, I  # X,V,I are returned after all iterations are complete


def main():
    n = 100  # spatial grid size.
    M = 5  # number of electrons injected per turn.
    nk = 500  # number of turns to simulate.
    u0 = 5  # threshold velocity.
    p = 0.25  # probability that ionization will occur
    Msig = 1

    if(len(sys.argv) > 1):
        n = int(sys.argv[1])
        M = int(sys.argv[2])
        nk = int(sys.argv[3])
        u0 = int(sys.argv[4])
        p = float(sys.argv[5])
        Msig = float(sys.argv[6])

    # Simulating the tubelight with given parameters
    X, V, I = Tubelight_sim(n, M, nk, u0, p, Msig)

    # Plotting Electron density Histogram
    plt.figure(0)
    plt.hist(X, bins=np.arange(1, n), ec='black')
    plt.ylabel(r'Electron density$\rightarrow$', fontsize=13)
    plt.xlabel(r'$x$$\rightarrow$', fontsize=13)
    plt.title("Electron Density Histogram")
    plt.show()

    # Plotting Light intensity histogram
    plt.figure(1)
    plt.hist(I, bins=np.arange(1, n), ec='black')
    plt.ylabel(r'$I\rightarrow$', fontsize=13)
    plt.xlabel(r'$x$$\rightarrow$', fontsize=13)
    plt.title("Light Intensity Histogram ")
    plt.show()

    # Plotting the electron phase space
    plt.figure(2)
    plt.plot(X, V, 'x')
    plt.title("Electron Phase Space")
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    plt.show()

    # Printing the table of xpos and count
    ints, bins, notused = plt.hist(I, bins=n)
    xpos = 0.5*(bins[0:-1]+bins[1:])
    print("Intensity Data")
    print("xpos \t count")
    for i in range(len(ints)):
        print(format(xpos[i], ".3f"), '\t', ints[i])


if __name__ == '__main__':
    main()
