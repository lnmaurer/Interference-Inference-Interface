# Interference Inference Interface
![The graphical user interface.](https://mywebspace.wisc.edu/lnmaurer/web/iii/double_slit.png)

This is a simple 2D FDTD simulation of an electromagnetic TMz wave made for teaching in an instructional lab.

### Interface

Consider the upper left plot in the interface; it shows Ez. Black is the most negative value and white is the most positive value. An plane wave comes in from the left of the barrier (the vertical red line). This wave has a fixed wavelength of 20 grid cells. Rather than simulate the incoming wave, it is calculated analytically. At the start of the simulation, the magnitude of this wave is ramped up gradually, to avoid potentially unstable high frequency components.

The barrier is a perfect electric conductor, and the openings in the barrier are hard sources that inject the waves in to the FDTD domain, the area to the right of the barrier. The other three sides of the FDTD domain are terminated with split-field perfectly matched layers. These reduce reflections from the boundaries to imperceptible levels, effectively giving the simulation open boundaries.

The bottom left plot shows Ezrms by default; it can also be set to show Ezrms^2 (proportional to intensity). In both cases, black is zero, and white is the maximum value. When the simulation starts or the barrier is modified, a timer appears over this plot, counting down until steady state is reached. After that, the averaging of Ezrms restarts to remove transients, and another timer appears, counting down until one time period is averaged.

Horizontal and vertical dashed lines are visible in both of these plots. Ez (yellow) and sqrt(2)*Ezrms (green, as an envelope) are shown along these vertical and horizontal slices in the remaining three plots (the two vertical plots are identical). The distances from the centers of the barrier's openings to the intersection of the slices is shown in cyan on the 2D plots and in the barrier control frame to the right.

x, y, Ez, and Ezrms at the point the slices intersect are shown in the center right area between the plots for the vertical slice. This is useful for homing in on maxima and minima; the slices can be moved either with the mouse or keyboard.

Finally, the openings in the barrier can be added, removed, and modified in the barrier control frame at the right of the interface.

The instructional lab computers (with Core 2 Duo processors) take ~55ms per simulation timestep, which results in a smooth animation. The simulation can be sped up using fast forward mode; it saves time by not updating the plots and runs until the current countdown is completed.

### Instructional Laboratory

For an example of a lab that uses this simulation, see the [worksheet](https://github.com/lnmaurer/Interference-Diffraction-Worksheet) I used.

### Code

The software is written in Python using NumPy for the calculations, TkInter for the interface, and the Python Imaging Library to produce the 2D plots. Those libraries are available for Windows, OS X, and Linux, so the program can run on any of those platforms. I will post executables presently.

### Notes
-  There's an executable for Window avaible one the [downloads page](https://github.com/lnmaurer/Interference-Inference-Interface/downloads)
  - I'll make other downloads available eventually
  - The executable makes a folder when it's run. It'll delete the folder when it exits if it exits normally. If it doesn't exit normally, feel free to delete the folder and the stuff in it.
-  This code works with Python 2.7
-  This program requires three libraries
  -  [Numpy](http://numpy.scipy.org/)
  -  TkInter, which is part of the default Python install
  -  [Python Imaging Library](http://www.pythonware.com/products/pil/). The only available downloads on that page are for Windows, but it as available on other platforms through other means. For example, on Mac OS X, it's [available](https://trac.macports.org/browser/trunk/dports/python/py-pil/Portfile) through MacPorts.

### Contact

If you have any questions, don't hesistate to contact me (@lnmaurer). For more information about me, also see my [webpage](https://mywebspace.wisc.edu/lnmaurer/web/)

### The code is available under the [GNU Public License Version 2](http://www.gnu.org/licenses/gpl-2.0.html)

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; using version 2
  of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.