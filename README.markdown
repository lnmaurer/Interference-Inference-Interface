# Interference Inference Interface
![The graphical user interface.](http://www.physics.wisc.edu/~lmaurer/code/iii/double_slit.png)
[Full-sized image](http://www.physics.wisc.edu/~lmaurer/code/iii/double_slit.png)

This is a simple 2D FDTD simulation of an electromagnetic TMz wave made for teaching in an instructional lab.

**[Here](http://www.youtube.com/watch?v=6QLj8_f9QYg)'s a short video introduction to the program.**

## Downloads

### Windows
[Here](http://www.physics.wisc.edu/~lmaurer/code/iii/iii.exe)'s a Window executable. Just double click it to run. Note that the executable makes a folder when it's run; it'll be deleted when if the program exits normally. If it doesn't exit normally, feel free to delete the folder and the stuff in it.

### OS X
[Here](http://www.physics.wisc.edu/~lmaurer/code/iii/iii.zip)'s an OS X executable in a zip file. Unzip it and then double click it to run.

Unfortunately, a known and apparently unfixable [bug](http://mail.python.org/pipermail/tkinter-discuss/2013-January/003343.html) in the interface library can make the interface unresponsive when the simulation is running on OS X. This code is using a workaround that seems to be effective, but it's possible you'll run in to problems. If so, don't hesitate to contact me.

Even with this workaround, the simulation runs more slowly -- making the animation more coppy -- on OS X.

### Linux
I've made a Debian package -- tested under Ubuntu 12.10 -- available [here](http://www.physics.wisc.edu/~lmaurer/code/iii/interference-inference-interface-69.deb).

If you're running something else, it's probably easy to get the needed libraries and get running. For example, Ubuntu 12.04 comes with python 2.7 by default, so ``sudo apt-get install python-numpy python-imaging python-imaging-tk`` will install the necessary libraries. After that, ```iii.py`` -- available [here](https://github.com/lnmaurer/Interference-Inference-Interface/zipball/master) -- can be run on the command line with ``python iii.py``

I understand a future version of Ubuntu will move to python 3, so you may need to add some '2.7's in there. E.g. ``sudo apt-get install python2.7 python2.7-numpy python2.7-imaging python2.7-imaging-tk`` then run it with ``python2.7 iii.py``.

## Interface

The interface shows a TMz wave with five plots. The two large plots show Ez (upper plot) and EzRMS (lower plot) throughout the simulation's domain. For both, black is the smallest value and white is the largest value, with shades of gray in between. The three smaller plots show Ez and +/-sqrt(2)*EzRMS -- an envelope for Ez -- along the horizontal and vertical dashed lines through the two larger plots. Those lines can be moved with the keyboard or mouse, and x, y, Ez, and EzRMS at those lines' intersection is displayed in the center right area between the two vertical plots. Knowing EzRMS at that point allows users to home in on extrema.

A plane wave -- with a wavelength of 20 grid cells -- enters from the left. Rather than simulate the incoming wave, it is calculated analytically. At the start of the simulation, the wave's magnitude is ramped up gradually to avoid potentially unstable high frequency components.

The barrier -- visible in both large plots -- is a perfect electric conductor, and the openings in the barrier are hard sources that inject the incoming wave in to the FDTD domain, the area to the right of the barrier. Split-field perfectly matched layers terminate the other three sides of the FDTD domain. This boundary reduces reflections to imperceptible levels, effectively giving the simulation open boundaries.

When the simulation starts or the barrier is modified, a timer appears over the EzRMS plot -- counting down until steady state is reached. Afterwards, EzRMS is reset to remove transients, and another timer appears for one time period of the wave. EzRMS is averaged over that time.

Openings in the barrier can be added, removed, and modified using the barrier control frame -- at the right of the interface.

Among its other features, the simulation has also has a fast forward mode, which saves time by not updating the plots. When in that mode, the simulation runs until the current countdown is done.

The laboratory's computers -- running Windows with Intel Core 2 Duo processors -- take about 55ms per simulation timestep, resulting in a smooth animation.

## Instructional Laboratory

For an example of a lab that uses this simulation, see the [worksheet](https://github.com/lnmaurer/Interference-Diffraction-Worksheet) I used.

## Code

The software is written in Python using NumPy for the calculations, TkInter for the interface, and the Python Imaging Library to produce the 2D plots. Those libraries are available for Windows, OS X, and Linux.

## Notes
-  This code works with Python 2.7; I'll probably make it compatible with Python 3 eventually -- once all the libraries have Python 3 versions.
-  This program requires three libraries
  -  [Numpy](http://numpy.scipy.org/)
  -  TkInter, which is part of the default Python install
  -  [Python Imaging Library](http://www.pythonware.com/products/pil/). The only available downloads on that page are for Windows, but it as available on other platforms through other means. For example, on Mac OS X, it's available through package managers like pip and MacPorts.

## Contact

If you have any questions, don't hesitate to contact me (@lnmaurer). For more information about me, see my [webpage](http://www.physics.wisc.edu/~lmaurer/)

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