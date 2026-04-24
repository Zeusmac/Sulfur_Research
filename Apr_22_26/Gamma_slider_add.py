import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# -------------------------
# Load ROOT file
# -------------------------
file = uproot.open("/Users/macwheeler/Sulfur_Research/HistRootFiles/dtime44S_main.root")
h2 = file["dtime"]

data = h2.values()   # shape: (Y, X) usually

# -------------------------
# Safety: ensure orientation
# -------------------------
# We assume:
# Y = time/slices
# X = energy channels

# -------------------------
# Parameters
# -------------------------
window = 10
initial_bin = 1001
n_slices = 500

# -------------------------
# Update function
# -------------------------
def update(val):

    i = int(slider.val)

    y = projections[i]

    line.set_ydata(y)

    # update x if needed (after cutting bins)
    x = np.arange(len(y))
    line.set_xdata(x)

    # dynamic autoscale
    ax.relim()
    ax.autoscale_view()

    ax.set_title(f"Slice {i}")

    fig.canvas.draw_idle()

cut_x = 20   # remove first 20 X bins

# -------------------------
# Build Y-projections (ProjectionY equivalent)
# -------------------------
projections = []

y2 = initial_bin

for i in range(n_slices):

    y1 = initial_bin
    y2 += window

    # sum over Y → produce X spectrum
    proj = np.sum(data[y1:y2, :], axis=0)

    # -------------------------
    # CUT FIRST 20 X BINS
    # -------------------------
    proj = proj[cut_x:]

    projections.append(proj)

projections = np.array(projections)

# -------------------------
# Plot setup
# -------------------------
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)

x = np.arange(projections.shape[1])

line, = ax.plot(x, projections[0])

ax.set_title("Y-Projection Spectra (cut first 20 bins)")
ax.set_xlabel("Channel (cut)")
ax.set_ylabel("Counts")

# -------------------------
# Slider
# -------------------------
ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, "Slice", 0, n_slices - 1, valinit=0, valstep=1)

slider.on_changed(update)

plt.show()
