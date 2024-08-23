import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey

Sankey(flows=[1, -0.1, -0.15, -0.05,-0.01,-0.12,-0.57],
   labels=['INLET', 'Stator', 'Rotor', 'Mech.', 'Th. loss', 'El.loss', 'El. Power OUT'],
   orientations=[ 0, 1, 1, 1,1, 1, 0]).finish()

# remove the frame in the sankey diagram:
a = plt.gca()
a.set_frame_on(False)

# Save it to pdf:
plt.savefig("sankey.pdf", bbox_inches='tight')