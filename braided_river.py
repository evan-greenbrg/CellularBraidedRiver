import copy
import math
import numpy as np
from matplotlib import pyplot as plt
import scipy.ndimage as ndimage


def range_check(i, j):
    if i == 0:
        ilow = i
        ihigh = i + 1
    elif i == HEIGHT -1:
        ilow = i - 1
        ihigh = i
    else:
        ilow = i - 1
        ihigh = i + 1

    if j == 0:
        jlow = j
        jhigh = j + 1
    elif j == WIDTH - 1:
        jlow = j - 1
        jhigh = j
    else:
        jlow = j - 1
        jhigh = j + 1

    return ilow, ihigh, jlow, jhigh


def get_upstream_slopes(elev, i, j):
    if i == 0:
        return np.zeros((3,))

    cell_elev = elev[i, j]
    slopes = np.zeros((3,))
    if j != 0:
        slopes[0] = (elev[i-1, j-1] - cell_elev) / (2**(1/2))
    if j != WIDTH - 1:
        slopes[2] = (elev[i-1, j+1] - cell_elev) / (2**(1/2))

    slopes[1] = (elev[i-1, j] - cell_elev) / (2**(1/2))

    return slopes


def get_upstream_water(inwater, i, j):
    if i == 0:
        return np.zeros((3,))

    upstream_water = np.zeros((3,))
    if j != 0:
        upstream_water[0] = inwater[i -1, j -1]
    if j != WIDTH - 1:
        upstream_water[2] = inwater[i - 1, j + 1]

    upstream_water[1] = inwater[i - 1, j]

    return upstream_water 


def get_downstream_slopes(elev, i, j):
    if i == HEIGHT - 1:
        return np.zeros((3,))

    cell_elev = elev[i, j]
    slopes = np.zeros((3,))
    if j != 0:
        slopes[0] = (cell_elev - elev[i+1, j-1]) / (2**(1/2))
    if j != WIDTH - 1:
        slopes[2] = (cell_elev - elev[i+1, j+1]) / (2**(1/2))

    slopes[1] = (cell_elev - elev[i+1, j]) / (2**(1/2))

    return slopes


def get_water_path(raw_slopes, cell_water):
    slopes = copy.deepcopy(raw_slopes)
    if (not len(slopes[slopes > 0])) and (len(slopes[slopes < 0])):
#        slopes[slopes == 0] = 1e-250
        where_zero = np.where(slopes == 0)[0]
        slopes = slopes[slopes != 0]
        water_path = np.multiply(
            cell_water,
            np.divide(
                np.abs(slopes)**-N, 
                np.sum(np.abs(slopes)**-N)
            )
        )
        water_path = np.insert(water_path, where_zero, 0)
    elif len(slopes[slopes > 0]) > 0:
        slopes[slopes < 0] = 0
        water_path = np.multiply(
            cell_water,
            np.divide(
                slopes**N, 
                np.sum(slopes**N)
            )
        )
    else:
        water_path = [cell_water / 3 for i in slopes]

    return water_path


def distribute(flux, grid, i, j):

    if j != 0:
        grid[j - 1] += flux[0]

    grid[j] += flux[1]

    if j != WIDTH - 1:
        grid[j + 1] += flux[2]

    return grid


def route_water_and_sediment(downstream_slopes, new_water_row,
                             cell_water, sediment, i, j):
    # Get outgoing water
    outgoing_water = get_water_path(downstream_slopes, cell_water)

    # Get outgoing sediment
    stream_power = (
        outgoing_water 
        * (np.add(downstream_slopes, CS))
    )
    stream_power[stream_power < 0] = 0
    outgoing_sediment = K * (stream_power**M)

    sediment[i, j, 0] = np.sum(outgoing_sediment)

    # Distribute water
    if i == HEIGHT - 1:
        return new_water_row, sediment 

    new_water_row = distribute(outgoing_water, new_water_row, i, j)

    # Distribute sediment
    sediment[i, :, 1] = distribute(
        outgoing_sediment, 
        sediment[i, :, 1], 
        i, j
    )

    return new_water_row, sediment


def main_loop(elev, water, sediment):
    for i, row in enumerate(elev):
        new_water_row = np.zeros((len(row)))

        cell_waters = []
        for j, val in enumerate(row):
            downstream_slopes = get_downstream_slopes(elev, i, j)
    #        upstream_slopes = get_upstream_slopes(elev, i, j)
    #        upstream_water = get_upstream_water(water, i, j)

            cell_water = water[i, j]
            new_water_row, sediment = route_water_and_sediment(
                downstream_slopes,
                new_water_row,
                cell_water,
                sediment,
                i, j
            )

            cell_waters.append(cell_water)
        row_water = round(np.sum(water[i,:]), 0)
        new_row_water = round(np.sum(new_water_row), 0)

        if i != HEIGHT - 1:
            water[i + 1, :] = new_water_row

    return water, sediment


# 1. Initialize Grid
NIT = 1000
E = .35
N = 0.5
CS = 300000
K = 10**-22
M = 2.5

WIDTH = 22
HEIGHT = 500
elev = np.zeros((HEIGHT, WIDTH))
water = np.zeros((HEIGHT, WIDTH)) # 0 - incoming 1 - outgoing
sediment = np.zeros((HEIGHT, WIDTH, 2)) # 0 - outgoing 1 -incoming 

# 2. Initialize slope
slope = 100000
e0 = 0 
initial_slope = [e0 for i in np.arange(1, 501)] + np.arange(1, 501) * slope
for i, column in enumerate(elev.T):
    elev[:, i] = initial_slope

# 3. Add noise
elev = np.flip(elev)
slope_elev = copy.deepcopy(elev)
elev_diff = elev[0, 0] - elev[1, 0]
noise = np.random.normal(0, scale=elev_diff, size=(HEIGHT, WIDTH))
elev = np.add(elev, noise)

# 4. set incoming water
q0 = 10000
input_loc = [int(WIDTH/2)]
q0 = q0 / len(input_loc)
for loc in input_loc:
    water[0, loc] = q0

# 5. Actually start model
# BLOWS UP PROBLEM!
for i in range(NIT):
    print(i)
    water, sediment = main_loop(elev, water, sediment)
    # Evolve the bed
    delev = np.subtract(sediment[:, :, 1], sediment[:, :, 0])
    print(np.max(delev))
    print(np.min(delev))
    elev += delev
    if math.isnan(elev[0,0]):
        break






test = np.sum(grid[:,:,1], axis=1)
plt.imshow(grid[:,:,1])
plt.show()


fig, axs = plt.subplots(1,2, sharex=True, sharey=True)
axs[0].imshow(water[:,:,0])
axs[1].imshow(water[:,:,1])
plt.show()

plt.imshow(water[:,:,1])
plt.show()
axs[0].imshow(water[:,:,0])


plt.imshow(water[:,:,0] == water[:,:,1])

test = np.sum(water[:,:,1], axis=1)
x = [i for i, val in enumerate(test)]
plt.plot(x, test)
plt.show()

