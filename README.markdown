# Interference Inference Interface
![The graphical user interface.](https://mywebspace.wisc.edu/lnmaurer/web/iii/double_slit.png)

This is a simple 2D FDTD simulation of an electromagnetic TMz wave made for teaching in an instructional lab.

[Full size picture of the interface](https://mywebspace.wisc.edu/lnmaurer/web/iii/double_slit.png)

## Downloads

### Windows
There's an executable for Window available on the [downloads page](https://github.com/lnmaurer/Interference-Inference-Interface/downloads). Just double click it to run. Note that The executable makes a folder when it's run; it'll be deleted when if the program exits normally. If it doesn't exit normally, feel free to delete the folder and the stuff in it.

### OS X
I'll try making an OS X version once [this bug](http://bugs.python.org/issue15574) in the OS X version of TkInter is fixed. But for now, the OS X version of TkInter is too buggy.

### Linux
I haven't made an executable for linux yet because it's easy to get the needed libraries. For exaple, Ubuntu 12.04 comes with python 2.7 by default, so ``sudo apt-get install python-numpy python-imaging`` will install the necessary libraries. After that, ```iii.py`` -- available [here](https://github.com/lnmaurer/Interference-Inference-Interface/zipball/master) -- can be run on the command line with ``python iii.py``

I understand the next version of Ubuntu will move to python 3, so you may need to add some '2.7's in there. E.g. ``sudo apt-get install python2.7 python2.7-numpy python2.7-imaging`` then run it with ``python2.7 iii.py``.

## Interface

The interface shows a TMz wave with five plots. The two large plots show Ez (upper plot) and EzRMS (lower plot) throughout the simulation's domain. For both, black is the smallest value and white is the largest value, with shades of gray in between. The three smaller plots show Ez and +/-sqrt(2)*EzRMS -- an envelope for Ez -- along the horizontal and vertical dashed lines through the two larger plots. Those lines can be moved with the keyboard or mouse, and x, y, Ez, and EzRMS at those lines' intersection is displayed in center right area between the two vertical plots. Knowing EzRMS at that point allows users to home in on extrema.

A plane wave -- with a wavelength of 20 grid cells -- enters from the left. Rather than simulate the incoming wave, it is calculated analytically. At the start of the simulation, the wave's magnitude is ramped up gradually to avoid potentially unstable high frequency components.

The barrier -- visible in both large plots -- is a perfect electric conductor, and the openings in the barrier are hard sources that inject the incoming wave in to the FDTD domain, the area to the right of the barrier. Split-field perfectly matched layers terminate the other three sides of the FDTD domain. This boundary reduces reflections to imperceptible levels, effectively giving the simulation open boundaries.

When the simulation starts or the barrier is modified, a timer appears over the EzRMS plot -- counting down until steady state is reached. Afterwards, EzRMS is reset to remove transients, and another timer appears for one time period of the wave. EzRMS is averaged over that time.

Openings in the barrier can be added, removed, and modified using the barrier control frame -- at the right of the interface.

Among its other features, the simulation has also has a fast forward mode, which saves time by not updating the plots. When in that mode, the simulation runs until the current countdown is done.

The laboratory's computers -- running Windows with Intel Core 2 Duo processors -- take about 55ms per simulation timestep, resulting in a smooth animation.

## Instructional Laboratory

For an example of a lab that uses this simulation, see the [worksheet](https://github.com/lnmaurer/Interference-Diffraction-Worksheet) I used.

## Code

The software is written in Python using NumPy for the calculations, TkInter for the interface, and the Python Imaging Library to produce the 2D plots. Those libraries are available for Windows, OS X, and Linux, so the program can run on any of those platforms. I will post executables presently.

## Notes
-  This code works with Python 2.7; I'll probably make it compatible with Python 3 eventually -- once all the libraries have Python 3 versions.
-  This program requires three libraries
  -  [Numpy](http://numpy.scipy.org/)
  -  TkInter, which is part of the default Python install
  -  [Python Imaging Library](http://www.pythonware.com/products/pil/). The only available downloads on that page are for Windows, but it as available on other platforms through other means. For example, on Mac OS X, it's [available](https://trac.macports.org/browser/trunk/dports/python/py-pil/Portfile) through MacPorts.

## Contact

If you have any questions, don't hesitate to contact me (@lnmaurer). For more information about me, see my [webpage](https://mywebspace.wisc.edu/lnmaurer/web/)

## The code is available under the [GNU Public License Version 2](http://www.gnu.org/licenses/gpl-2.0.html)

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