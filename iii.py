import Tkinter, tkFileDialog, ttk, tkMessageBox
import csv #for exporting in CSV
import time #for testing how long steps take
from numpy import * #so that we can use numpy.round, which otherwise conflicts with python's built in rounding function
import numpy
from PIL import Image, ImageTk, ImageDraw

#CONSTANTS---------------------------------------------------------------------
mu0      = 1.2566370614e-6    #Vacuum permeability
epsilon0 = 8.854187817620e-12 # Vacuum permittivity
c        = 1/(mu0*epsilon0)**0.5
  
Nx = 600 #width of simulation and view canvas
Ny = 300 #height of simulation and view canvas
barrierX = 100 #x position of the barrier
  
plotD = 100 #dimension of plot
  
d  = 2.0/Nx #spatial grid element size
dt = d/c/2**0.5 #time step -- this choice is good for a 2D simulation in vacuum

#the wave:
lamb  = 20*d		#wavelength -- 20 points per wavelength yeilds good results
k     = 2*pi/lamb	#wavenumber
omega = 2*pi*c/lamb	#angular frequency
tau   = lamb/c		#time period

#stability time step counts
nStable = int((barrierX + sqrt((Nx-barrierX)**2 + Ny**2))*d/c/dt) #time steps until Ez beocmes stable: to barrier then longest way across right domain
nAvgStable = int(10*tau/dt) #time steps of averaging to get the average stable (once Ez is already stable)

#Yee algorithm update coefficients
Db = dt/mu0/d
Cb = dt/epsilon0/d

#Mur RBC update coefficients
Ma = (c*dt - d)/(c*dt + d)
Mb = 2*d/(c*dt + d)
Mc = (c*dt)**2/2/d/(c*dt + d)
  
NEZx = Nx #number of Ez grid points in x direction -- goes all the way to the edge
NEZy = Ny #number of Ez grid points in y direction -- goes all the way to the edge
#range updated by Yee algorithm:
EZx_range = slice(barrierX+1, NEZx-1) #right end is RBC and left side is explicitly calculated, so not updated using Yee algorithm
EZy_range = slice(1, NEZy-1) #top and bottom are RBCs, so not updated using Yee algorithm
#range updated by explicit formula:
EZx_ex_range = slice(0, barrierX+1)
EZy_ex_range = slice(0, NEZy)
#everything
EZx_all = slice(0,NEZx)
EZy_all = slice(0,NEZy)
  
NHXx = Nx     #Hx goes all the way to x endpoints
NHXy = Ny - 1 #Hx not at y endpoints, so -1
HXx_range = slice(0, NHXx) #all positions updated using Yee
HXy_range = slice(0, NHXy)

NHYx = Nx - 1 #Hy not at x endpoints, so -1
NHYy = Ny     #Hy goes all the way to y endpoints
HYx_range = slice(0, NHYx) #all y positions updated using Yee
HYy_range = slice(0, NHYy) 

#THE METHODS-------------------------------------------------------------------
def exportData():
  fileName = tkFileDialog.asksaveasfilename(filetypes=[('CSV','*.csv')], title="Export data as...")
  if fileName != '': #'' is returned if the user hits cancel
    writer = csv.writer(open(fileName, "w"))
    writer.writerow(('x_slice',sliceX))
    writer.writerow(('y_slice',sliceY))
    writer.writerow(('time steps',n))
    writer.writerow(('time steps averaged',nAveraging))
    
    writer.writerow(('Ez'))
    row = ['y\\x']
    row.extend(range(0,Nx))
    writer.writerow(row)
    for i, r in enumerate(Ez[:,:,0]):
      row = [i]
      row.extend(r)
      writer.writerow(row)
      
    writer.writerow(('Ez_RMS'))
    row = ['y\\x']
    row.extend(range(0,Nx))
    writer.writerow(row)
    for i, r in enumerate(EzRMS[:,:]):
      row = [i]
      row.extend(r)
      writer.writerow(row)
  
def addOpening():
  global gaps
  global nEnd
  
  bot = int(bottomEntry.get())
  top = int(topEntry.get())
  if bot > 0 and top < Ny and bot < top: #some initial checks that top and bot better pass
    newOrder = gaps[:]
    newOrder.append([bot,top])
    newOrder.sort()
    newOrderFlat = [ypos for pair in newOrder for ypos in pair] #flatten
    if (newOrderFlat == sorted(newOrderFlat) #catches things like [50,100],[90,120] -- overlapping openings
      and len(list(set(newOrderFlat))) == len(newOrderFlat)): #catchs things like [50,100],[100,120] -- openings with no space between them
      gaps = newOrder
      nEnd = n + nStable
      redrawBarrierFrame()
      conditionalRedraw()
  
def redrawBarrierFrame():
  global barrierFrames
  global strVars
  global distStrVars
  
  gaps.sort()
  for f in barrierFrames:
    f.destroy() #get rid of old frames
  barrierFrames = []
  strVars       = [] #need to save StringVars or else they get garbage collected
  distStrVars   = [] #holds stringVars that report the distance from the center of the gaps to the slice intersection
  r = 3 #rows 0,1,2 already taken by widgets for adding an opeing
    
  for gap in gaps:
    bn = r-3 #barrier number
    frame = ttk.Labelframe(barrierFrame, text="Opening {}".format(bn))
    frame.grid(column=0, row=r, columnspan=2, sticky='nesw', padx=5, pady=5)
      
    top = Tkinter.IntVar()
    bottom = Tkinter.IntVar()
    distStrVar = Tkinter.StringVar()
    ttk.Label(frame, text="Bottom:").grid(column=0, row=0, sticky='nes', padx=5, pady=5)
    #having the following work is kind of tricky; the default parameter in the lambda is critical. See <http://mail.python.org/pipermail/tutor/2005-November/043360.html>
    entry = Tkinter.Spinbox(frame, width=4, textvariable=bottom, from_=0, to=Nx, command=lambda n=bn, tv=bottom: updateBarrierBottom(n,tv))
    entry.bind("<Return>",lambda arg, n=bn, tv=bottom: updateBarrierBottom(n,tv))
    entry.grid(column=1, row=0, sticky='nsw', padx=5, pady=5)
    ttk.Label(frame, text="Top:").grid(column=0, row=1, sticky='nes', padx=5, pady=5)
    entry = Tkinter.Spinbox(frame, width=4, textvariable=top, from_=0, to=Nx, command=lambda n=bn, tv=top: updateBarrierTop(n,tv))
    entry.bind("<Return>",lambda arg, n=bn, tv=top: updateBarrierTop(n,tv))
    entry.grid(column=1, row=1, sticky='nsw', padx=5, pady=5)
    ttk.Label(frame, textvariable=distStrVar).grid(column=0, row=2, sticky='nes', padx=5, pady=5)
    ttk.Button(frame, text='Remove', command=lambda n=bn: removeBarrier(n)).grid(column=0, row=3, sticky='nsew', columnspan=2, padx=5, pady=5)
    top.set(gap[1])
    bottom.set(gap[0])
    strVars.append(top)
    strVars.append(bottom) 
    distStrVars.append(distStrVar)
      
    barrierFrames.append(frame)
    r = r + 1
  updateDistStrVars()
      
def updateBarrierTop(barrierNumber, intVar):
  global gaps
  global nEnd
  
  value = intVar.get()
  if ((barrierNumber == (len(gaps)-1)) and (value < Ny) and (value > gaps[-1][0])) or ((barrierNumber < (len(gaps)-1)) and (value < gaps[barrierNumber+1][0]) and (value > gaps[barrierNumber][0])):
    gaps[barrierNumber][1] = value
    nEnd = n + nStable
    conditionalRedraw()
  else:
    intVar.set(gaps[barrierNumber][1])

def updateBarrierBottom(barrierNumber, intVar):
  global gaps
  global nEnd

  value = intVar.get()
  if ((barrierNumber == 0) and (value > 0) and (value < gaps[0][1])) or ((barrierNumber > 0) and (value > gaps[barrierNumber-1][1]) and (value < gaps[barrierNumber][1])):
    gaps[barrierNumber][0] = value
    nEnd = n + nStable
    conditionalRedraw()
  else:
    intVar.set(gaps[barrierNumber][0])    
  
def updateDistStrVars():
  for gap, tv in zip(gaps, distStrVars):
    tv.set("dist=" + str(int(round(sqrt(((gap[0]+gap[1])/2.0-sliceY)**2+(barrierX-sliceX)**2)))) + "d")
    
def invertedGaps():
  """Returns [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]"""
  invGaps = [ypos for pair in gaps for ypos in pair] #flatten gaps
  invGaps.insert(0,0) #insert zero at the 0th position
  invGaps.append(Ny) #put Ny for the last item
  return invGaps
  
def removeBarrier(barrierNumber):
  global gaps
  global nEnd

  del gaps[barrierNumber]
  nEnd = n + nStable
  redrawBarrierFrame()
  conditionalRedraw()
    
def clearCanvasBindings(eventObj):
  Ezcanvas.bind("<B1-Motion>", lambda e: None)
  HorizPlotCanvas.bind("<B1-Motion>", lambda e: None)
  EzRMScanvas.bind("<B1-Motion>", lambda e: None)
  
def updateEzPlot():
  global ezPlot
  #notes:
  #1)Image.fromstring can only handle 32bit floats, so need to do that conversion
  #2)0 (and below) are black, 255 and above are white, shades of gray inbetween
  #3)need to store array in (height, width) format
  data = float32((transpose(Ez[:,:,0])/maxY + 1)/2*256) #+1 so that zero is in the center
  im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
  ezPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas
    
def updateEzRMSPlot():
  global ezRMSPlot
  
  if avgSetting.get() == 'sq': #want to display Ez_RMS^2
    data = 256*(float32(transpose(EzSQ)/nAveraging/maxRMSY**2))  
  else: #want to display Ez_RMS
    data = 256*(float32(transpose(EzRMS)/maxRMSY))    
  im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
  ezRMSPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas

def makeAvgedTrace(xS, yS, invert=False, othercoord=None):
  if avgSetting.get() == 'amp': #want to display amplitude = EzRMS*sqrt(2) with zero in middle of plot
    v = int_(numpy.round((1+EzRMS[xS,yS]*sqrt(2)/maxY)*(plotD/2)))
    if othercoord != None: #in this case, let's plot +/-amplitude, not just amplitude
      othercoord.extend(othercoord[::-1]) #extend the coordinates so they're like [0,1...,Nx-1,Nx,Nx,Nx-1,...,1,0]
      v = concatenate((v,100-v[::-1])) #extend the values so they're like [V0,...,Vn,-Vn,...,-V0]
  elif avgSetting.get() == 'rms': #want to display EzRMS with zero at bottom of plot
    v = int_(numpy.round((EzRMS[xS,yS]/maxRMSY)*(plotD-1)))
  else: #want to display EzRMS^2 with zero at bottom of plot
    v = int_(numpy.round((EzSQ[xS,yS]/nAveraging/maxRMSY**2)*(plotD-1)))
    
  if invert: #mirror results across axis
    v = plotD - v
  return v
    
def plot(draw, x, y, color):
  """Plots the data on the drawing in the given color"""
  if traceSetting.get() == 'line': #line plot
    draw.line(zip(x,y), fill=color)
  else: #dot plot
    draw.point(zip(x,y), fill=color)
      
def updateHorizPlot():
  global horizPlot
  
  im = Image.new('RGB', (Nx,plotD))
  draw = ImageDraw.Draw(im)
  x = range(0,Nx)
  #plot EzRMS
  y = makeAvgedTrace(EZx_all, sliceY, invert=True, othercoord=x)
  plot(draw, x, y, 'green')
  #plot Ez
  y = int_(numpy.round((-Ez[:,sliceY,0]/maxY+1)*(plotD/2)))
  plot(draw, x, y, 'yellow')
  #turn plot in to a format the canvas can use
  horizPlot = ImageTk.PhotoImage(image=im)

def updateVertPlot():
  global vertPlot
  
  im = Image.new('RGB', (plotD, Ny))
  draw = ImageDraw.Draw(im)
  y = range(0,Ny)
  #plot EzRMS
  x = makeAvgedTrace(sliceX, EZy_all, othercoord=y)
  plot(draw, x, y, 'green')
  #plot Ez
  x = int_(numpy.round((Ez[sliceX,:,0]/maxY+1)*(plotD/2)))
  plot(draw, x, y, 'yellow')
  #turn plot in to a format the canvas can use  
  vertPlot = ImageTk.PhotoImage(image=im)  
    
def redrawCanvases():
  global maxRMSY
  global maxY
  global EzRMS
    
  #first, clear everything off the canvases (but don't delete the canvases themselves)
  Ezcanvas.delete('all')
  HorizPlotCanvas.delete('all')
  EzRMScanvas.delete('all')
  VertPlotCanvas1.delete('all')
  VertPlotCanvas2.delete('all')
    
  EzRMS = sqrt(EzSQ/(nAveraging+1e-24))
  maxRMSY = numpy.max(EzRMS)
  if maxRMSY == 0:
    maxRMSY = 1

  tempMaxY = numpy.max(abs(Ez))
  if tempMaxY > maxY:
    maxY = tempMaxY
  
  #now, put the plots on the canvases
  updateEzPlot()
  Ezcanvas.create_image(0,0,image=ezPlot,anchor=Tkinter.NW)
  updateHorizPlot()
  HorizPlotCanvas.create_image(0,0,image=horizPlot,anchor=Tkinter.NW)  
  updateEzRMSPlot()
  EzRMScanvas.create_image(0,0,image=ezRMSPlot,anchor=Tkinter.NW)
  HorizPlotCanvas.create_image(0,0,image=horizPlot,anchor=Tkinter.NW)  
  updateVertPlot()
  VertPlotCanvas1.create_image(0,0,image=vertPlot,anchor=Tkinter.NW)  
  VertPlotCanvas2.create_image(0,0,image=vertPlot,anchor=Tkinter.NW)  
  
  #next, draw the barrier
  invGaps = invertedGaps()
  for i in range(0,len(invGaps),2):
    Ezcanvas.create_line([(barrierX,invGaps[i]),(barrierX,invGaps[i+1])], width=1, fill='red')    
    EzRMScanvas.create_line([(barrierX,invGaps[i]),(barrierX,invGaps[i+1])], width=1, fill='red')
  
  #draw lines from the center of the gaps to the slice intersection, if desired
  if distLineSetting.get() == "lines":
    for gap, distStrVar in zip(gaps,distStrVars):
      #draw the lines
      yGap = (gap[0]+gap[1])/2
      Ezcanvas.create_line([(barrierX,yGap),(sliceX,sliceY)], width=1, fill='orange', dash='.')
      EzRMScanvas.create_line([(barrierX,yGap),(sliceX,sliceY)], width=1, fill='orange', dash='.')
      #draw text showing the distances
      yText = (yGap + sliceY)/2
      xText = (barrierX+sliceX)/2
      Ezcanvas.create_text(xText, yText, text=distStrVar.get(), fill="Cyan", font=("system", "12"))
      EzRMScanvas.create_text(xText, yText, text=distStrVar.get(), fill="Cyan", font=("system", "12"))
  
  #now, draw the horizontal slice
  lineID = Ezcanvas.create_line([(0,sliceY),(Nx,sliceY)], width=1, fill='yellow', dash='-') 
  Ezcanvas.tag_bind(lineID, "<Button-1>",  horizClickMethod)
  lineID = EzRMScanvas.create_line([(0,sliceY),(Nx,sliceY)], width=1, fill='green', dash='-')
  EzRMScanvas.tag_bind(lineID, "<Button-1>",  horizClickMethod)    
  lineID = VertPlotCanvas1.create_line([(0,sliceY),(plotD,sliceY)], width=1, fill='blue', dash='-')
  VertPlotCanvas1.tag_bind(lineID, "<Button-1>",  horizClickMethod)
  lineID = VertPlotCanvas2.create_line([(0,sliceY),(plotD,sliceY)], width=1, fill='blue', dash='-')
  VertPlotCanvas2.tag_bind(lineID, "<Button-1>",  horizClickMethod)    
    
  #the vertical slice
  lineID = Ezcanvas.create_line([(sliceX,0),(sliceX,Ny)], width=1, fill='yellow', dash='-') 
  Ezcanvas.tag_bind(lineID, "<Button-1>",  vertClickMethod)
  lineID = EzRMScanvas.create_line([(sliceX,0),(sliceX,Ny)], width=1, fill='green', dash='-')
  EzRMScanvas.tag_bind(lineID, "<Button-1>",  vertClickMethod)
  lineID = HorizPlotCanvas.create_line([(sliceX,0),(sliceX,plotD)], width=1, fill='blue', dash='-')
  HorizPlotCanvas.tag_bind(lineID, "<Button-1>",  vertClickMethod)
 
  #if we haven't reached steady state yet, display a countdown with the number of steps until stability on EzRMScanvas
  #this is drawn last so that it's on top of everything
  if n <= nEnd:
    EzRMScanvas.create_text(Nx/2,Ny/2, text=str(nEnd-n), fill="Red", font=("system", "75"))
  elif nAveraging < nAvgStable:
    EzRMScanvas.create_text(Nx/2,Ny/2, text=str(nAvgStable-nAveraging), fill="Magenta", font=("system", "75"))
 
  #finially, show the current information about the slice intersection point
  tStringVar.set("t=" + str(n) + "dt")
  EzStringVar.set("Ez={:+.4f}".format(Ez[sliceX,sliceY,0]))
  EzRMSStringVar.set("EzRMS={:.4f}".format(EzRMS[sliceX,sliceY]))
    
def horizClickMethod(eventObj):
  Ezcanvas.bind('<B1-Motion>', horizDragMethod)
  EzRMScanvas.bind('<B1-Motion>', horizDragMethod)
  VertPlotCanvas1.bind('<B1-Motion>', horizDragMethod)
  VertPlotCanvas2.bind('<B1-Motion>', horizDragMethod)
    
def horizDragMethod(eventObj):
  setSliceY(eventObj.y)

def conditionalRedraw():
  """In many cases, we want to redraw the canvas because something has changed. However, if it's running, the canvas will be redrawn soon anyway, so we don't need to do an extra redraw."""
  if not running:
    redrawCanvases()
      
def setSliceY(y):
  global sliceY
  
  if (y >= 0) and (y < Ny):
    sliceY = y
    yStringVar.set("y=" + str(sliceY) + "d")
    updateDistStrVars()
    conditionalRedraw()
    
def vertClickMethod(eventObj):
  Ezcanvas.bind('<B1-Motion>', vertDragMethod)
  EzRMScanvas.bind('<B1-Motion>', vertDragMethod)
  HorizPlotCanvas.bind('<B1-Motion>', vertDragMethod)
    
def vertDragMethod(eventObj):
  setSliceX(eventObj.x)
    
def setSliceX(x):
  global sliceX
  
  if (x >= 0) and (x < Nx):
    sliceX = x
    xStringVar.set("x=" + str(sliceX) + "d")
    updateDistStrVars()
    conditionalRedraw()  
    
def resetAveraging():
  global nAveraging
  global EzSQ
  
  nAveraging = 0
  EzSQ = zeros((NEZx, NEZy))
  conditionalRedraw()
      
def reset():
  global n
  global Ez
  global Hx
  global Hy
  global t
  global maxY
  global nEnd
  global running
  
  running = False
  Hx = zeros((NHXx, NHXy))
  Hy = zeros((NHYx, NHYy))
  n = 0
  nEnd = nStable
  Ez = zeros((NEZx, NEZy, 3))
  maxY = 1.0
  resetAveraging() #contains a conditionalRedraw()
  
def start():
  global running
  
  if not running:
    running = True
    run()
    
def stop():
  global running

  running = False
  
def fastForward():
  global nEnd
  global fastForwarding
  global fastForwardWithAvg
  
  fastForwarding = True
  if nEnd < n: #we've already reached stability, so now we're fast forwarding to get good averaging
    fastForwardWithAvg = True
  else:
    fastForwardWithAvg = False
  root.after(1,fastForwardStep)
    
def fastForwardStep():
  global ezPlot
  global fastForwarding

  if (fastForwardWithAvg and nAveraging < nAvgStable) or n <= nEnd:
    #clear canvases and display countdown
    Ezcanvas.delete('all')
    EzRMScanvas.delete('all')
    if fastForwardWithAvg:
      Ezcanvas.create_text(Nx/2,Ny/2, text=str(nAvgStable-nAveraging), fill="Magenta", font=("system", "75")) 
      EzRMScanvas.create_text(Nx/2,Ny/2, text=str(nAvgStable-nAveraging), fill="Magenta", font=("system", "75"))       
    else:
      Ezcanvas.create_text(Nx/2,Ny/2, text=str(nEnd-n), fill="Red", font=("system", "75")) 
      EzRMScanvas.create_text(Nx/2,Ny/2, text=str(nEnd-n), fill="Red", font=("system", "75")) 
    step(avg=fastForwardWithAvg) #only average if we're not going to reset right after we're done fast forwarding
    root.after(1,fastForwardStep)
  else: #stop fast forwarding
    fastForwarding = False
    if not fastForwardWithAvg:
      resetAveraging()
    if running:
      run()
    else:
      start()
  
def run():
  if running:
    timer = time.clock()
    if n == nEnd:
      resetAveraging()
    step()
    redrawCanvases()
    print str(time.clock()-timer)
    if not fastForwarding:
      root.after(1,run)
  
def step(avg=True):
  global Ez
  global Hz
  global Hy
  global n
  global nAveraging
  global EzSQ
  
  #take care of Hx and Hy using the standard Yee algorithm
  Hx[HXx_range, HXy_range] = Hx[HXx_range, HXy_range] + Db*(Ez[HXx_range, HXy_range, 0] - Ez[HXx_range, 1:(NHXy+1), 0]) #HXy_range+1
  Hy[HYx_range, HYy_range] = Hy[HYx_range, HYy_range] + Db*(Ez[1:(NHYx+1), HYy_range, 0] - Ez[HYx_range, HYy_range, 0]) #HYx_range+1
    
  #move old Ez back in the array
  Ez = roll(Ez,1,axis=2)
    
  #do the normal Yee updates on Ez in the relevant range
  Ez[EZx_range, EZy_range, 0] = Ez[EZx_range, EZy_range, 1] + Cb*(Hy[EZx_range, EZy_range] - Hy[barrierX:(NEZx-2),EZy_range] + Hx[EZx_range, 0:NEZy-2] - Hx[EZx_range, EZy_range])
    
  #now take care of the three edge and two corner Mur RBCs
  #for x=NEZx-1
  rng   = slice(1,NEZy-1) #range of everything in y except the corners
  rngp1 = slice(2,NEZy)
  rngm1 = slice(0,NEZy-2)
  Ez[-1,rng,0] = -Ez[-2,rng,2] + Ma*(Ez[-2,rng,0] + Ez[-1,rng,2]) + Mb*(Ez[-1,rng,1] + Ez[-2,rng,1]) + Mc*(Ez[-1,rngp1,1] - 2*Ez[-1,rng,1] + Ez[-1,rngm1,1] + Ez[-2,rngp1,1] - 2*Ez[-2,rng,1] + Ez[-2,rngm1,1])

  #for y=0
  rng = slice(barrierX+2,NEZx-1) #range of everything in x except right corner
  rngp1 = slice(barrierX+3,NEZx)
  rngm1 = slice(barrierX+1,NEZx-2)
  Ez[rng,0,0] = -Ez[rng,1,2] + Ma*(Ez[rng,1,0] + Ez[rng,0,2]) + Mb*(Ez[rng,0,1] + Ez[rng,1,1]) + Mc*(Ez[rngp1,0,1] - 2*Ez[rng,0,1] + Ez[rngm1,0,1] + Ez[rngp1,1,1] - 2*Ez[rng,1,1] + Ez[rngm1,1,1])

  #for y=NEZy
  Ez[rng,-1,0] = -Ez[rng,-2,2] + Ma*(Ez[rng,-2,0] + Ez[rng,-1,2]) + Mb*(Ez[rng,-1,1] + Ez[rng,-2,1]) + Mc*(Ez[rngp1,-1,1] - 2*Ez[rng,-1,1] + Ez[rngm1,-1,1] + Ez[rngp1,-2,1] - 2*Ez[rng,-2,1] + Ez[rngm1,-2,1])

  #now for the corners
  Ez[-1,0,0] = Ez[-2,1,2] #bottom right
  Ez[-1,-1,0] = Ez[-2,-2,2] #top right
    
  #finially, update the time and the sources
  n += 1

  #next, update Ez in the calculated area -- x=0[0,barrierX]
  #x and y for points on the barrier and to the left
  x, y = d * mgrid[EZx_ex_range, EZy_ex_range]
  
  #want wave to ramp up the wave's magnitude gradually and propigate at speed of light
  #logistic growth makes it come in gradually and use of retarded time there and in step function enforces propigation
  #to increase efficency, don't worry about either after t > (barrierX*d/c + 10*tau)
  #at that point, the wave will already have reached the barrier and grown to nearly full strength
  Ez[EZx_ex_range, EZy_ex_range, 0] = sin(k*x-omega*n*dt)
  if n*dt < (barrierX*d/c + 10*tau):
    tr = n*dt - x/c
    Ez[EZx_ex_range, EZy_ex_range, 0] *= ((tr > 0).astype(float))/(1+exp(-(tr-3*tau)/tau))
  
  #now, enforce Ez=0 on barrier
  invGaps = invertedGaps()
  for i in range(0,len(invGaps),2):
    Ez[barrierX, invGaps[i]:invGaps[i+1], 0] = 0     
    
  #first, only average if avg is true
  #2nd, stop averaging once we've got enough points (i.e. after nAveraging == nAvgStable)
  #however, if Ez still isn't stable (i.e. n < nEnd) then keep averaging anyway
  if avg and (nAveraging < nAvgStable or n < nEnd):
    nAveraging += 1
    EzSQ += square(Ez[:,:,0])

#GLOBAL VARIABLES--------------------------------------------------------------
n = 0 #the current time step
nAveraging = 0 #number of time steps average has been running
running = False #the simulation isn't currently running
Ez = zeros((NEZx, NEZy,3)) #3rd dimension to keep track of two past values of Ez
EzSQ = zeros((NEZx, NEZy))
Hx = zeros((NHXx, NHXy))
Hy = zeros((NHYx, NHYy))
maxY = 1 #contains the largest Ez seen
fastForwarding = False
nEnd = nStable #number of steps stability will be reached at

#THE GUI-----------------------------------------------------------------------
#The root window
root = Tkinter.Tk()
root.title("Interference Inference Interface")

#The menubar and menus
menubar = Tkinter.Menu(root)

#the file menu
filemenu = Tkinter.Menu(menubar, tearoff=0)
#filemenu.add_command(label="Save Experiment As", accelerator="Ctrl+S")
#filemenu.add_command(label="Load Experiment", accelerator="Ctrl+O")
#filemenu.add_separator()
filemenu.add_command(label="Export Data", accelerator="Ctrl+E", command=exportData)
filemenu.add_separator()
filemenu.add_command(label="Exit", accelerator="Ctrl+Q", command=root.quit)
#bind keys to the actions
root.bind_all('<Control-e>', lambda arg: exportData())
root.bind_all('<Control-q>', lambda arg: root.quit())
menubar.add_cascade(label="File", menu=filemenu)

#the edit menu
editmenu = Tkinter.Menu(menubar, tearoff=0)
editmenu.add_command(label="Cut", accelerator="Ctrl+X", command=lambda: root.event_generate('<Control-x>'))
editmenu.add_command(label="Copy", accelerator="Ctrl+C", command=lambda: root.event_generate('<Control-c>'))
editmenu.add_command(label="Paste", accelerator="Ctrl+V", command=lambda: root.event_generate('<Control-v>'))
menubar.add_cascade(label="Edit", menu=editmenu)

#the view menu
viewmenu = Tkinter.Menu(menubar, tearoff=0)
traceSetting = Tkinter.StringVar(value="line")
viewmenu.add_radiobutton(label="Trace Setting:", state="disabled")
viewmenu.add_radiobutton(label="Lines", variable=traceSetting, value="line", command=conditionalRedraw)
viewmenu.add_radiobutton(label="Dots", variable=traceSetting, value="dot", command=conditionalRedraw)

viewmenu.add_separator()
avgSetting = Tkinter.StringVar(value='amp')
viewmenu.add_radiobutton(label="Displayed Average:", state="disabled")
viewmenu.add_radiobutton(label="Amplitude", variable=avgSetting, value="amp", command=conditionalRedraw)
viewmenu.add_radiobutton(label="Ez_rms", variable=avgSetting, value="rms", command=conditionalRedraw)
viewmenu.add_radiobutton(label="Ez_rms^2", variable=avgSetting, value="sq", command=conditionalRedraw)

viewmenu.add_separator()
distLineSetting = Tkinter.StringVar(value="lines")
viewmenu.add_checkbutton(label="Gap to slice lines", variable=distLineSetting, onvalue="lines", offvalue="nolines", command=conditionalRedraw)
menubar.add_cascade(label="View", menu=viewmenu)

#the simulation menu
simmenu = Tkinter.Menu(menubar, tearoff=0)
simmenu.add_command(label="Run", accelerator="Ctrl+R", command=start)
root.bind("<Control-r>", lambda arg: start())
simmenu.add_command(label="Stop", accelerator="Ctrl+S", command=stop)
root.bind("<Control-s>", lambda arg: stop())
simmenu.add_command(label="Reset Simulation", accelerator="Ctrl+Shift+R", command=reset)
root.bind("<Control-R>", lambda arg: reset())
simmenu.add_command(label="Reset Average", accelerator="Ctrl+Shift+A", command=resetAveraging)
root.bind("<Shift-A>", lambda arg: resetAveraging())
simmenu.add_command(label="Fast Forward", accelerator="Ctrl+F", command=fastForward)
root.bind("<Control-f>", lambda arg: fastForward())
menubar.add_cascade(label="Simulation", menu=simmenu)

#the help menu
helpmenu = Tkinter.Menu(menubar, tearoff=0)
helpmenu.add_command(label="Ask if you need help", state="disabled")
helpmenu.add_command(label="About", command=lambda: tkMessageBox.showinfo("About", "Leon's New Fangled Interference & Diffraction Simulator\n\nCopyright 2012 by Leon Maurer\n\nCode available under GNU Public License Version 2"))
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)

#The view frame
viewFrame = ttk.Labelframe(root, text='View')
viewFrame.grid(column=0, row=0, sticky='nsew',padx=5,pady=5)

Ezcanvas = Tkinter.Canvas(viewFrame, width=Nx, height=Ny)
Ezcanvas.grid(column=0, row=0, columnspan=3, rowspan=3, sticky='nsew', padx=5, pady=5)    

HorizPlotCanvas = Tkinter.Canvas(viewFrame, width=Nx, height=plotD)
HorizPlotCanvas.grid(column=0, row=3, columnspan=3, rowspan=5, sticky='nsew', padx=5, pady=5)

sliceY = 60 #position of horizontal slice
sliceX = Nx/2

xStringVar = Tkinter.StringVar(value="x=" + str(sliceX) + "d")
yStringVar = Tkinter.StringVar(value="y=" + str(sliceY) + "d")
tStringVar = Tkinter.StringVar()
EzStringVar = Tkinter.StringVar()
EzRMSStringVar = Tkinter.StringVar()
ttk.Label(viewFrame, textvariable=xStringVar).grid(column=3, row=3, sticky='nsew', padx=5, pady=0)
ttk.Label(viewFrame, textvariable=yStringVar).grid(column=3, row=4, sticky='nsew', padx=5, pady=0)
ttk.Label(viewFrame, textvariable=tStringVar).grid(column=3, row=5, sticky='nsew', padx=5, pady=0)    
ttk.Label(viewFrame, textvariable=EzStringVar).grid(column=3, row=6, sticky='nsew', padx=5, pady=0)    
ttk.Label(viewFrame, textvariable=EzRMSStringVar).grid(column=3, row=7, sticky='nsew', padx=5, pady=0)    

VertPlotCanvas1 = Tkinter.Canvas(viewFrame, width=plotD, height=Ny)
VertPlotCanvas1.grid(column=3, row=0, columnspan=1, rowspan=3, sticky='nsew', padx=5, pady=5)
VertPlotCanvas2 = Tkinter.Canvas(viewFrame, width=plotD, height=Ny)
VertPlotCanvas2.grid(column=3, row=8, columnspan=1, rowspan=3, sticky='nsew', padx=5, pady=5)

EzRMScanvas = Tkinter.Canvas(viewFrame, width=Nx, height=Ny)
EzRMScanvas.grid(column=0, row=8, columnspan=3, rowspan=3, sticky='nsew', padx=5, pady=5)

#make it so that, after dragging an element has ceased, the binding is reset so that further dragging won't move the element unless it gets clicked again first
#we do this by binding to any motion on any canvas
Ezcanvas.bind("<Motion>", clearCanvasBindings)
HorizPlotCanvas.bind("<Motion>", clearCanvasBindings)
EzRMScanvas.bind("<Motion>", clearCanvasBindings)

#make the arrow keys control the slice positions
root.bind("<Up>", lambda arg: setSliceY(sliceY-1))
root.bind("<Down>", lambda arg: setSliceY(sliceY+1))
root.bind("<Right>", lambda arg: setSliceX(sliceX+1))
root.bind("<Left>", lambda arg: setSliceX(sliceX-1))

#The barrier frame and intial barrier setup
barrierFrame = ttk.Labelframe(root, text='Barrier')
barrierFrame.grid(column=1,row=0,sticky='nsew',padx=5,pady=5)

ttk.Button(barrierFrame, text='Add Opening', command=addOpening).grid(column=0, row=0, columnspan=2, sticky='nsew', padx=5, pady=5)
ttk.Label(barrierFrame, text="Bottom:").grid(column=0, row=1, sticky='nes', padx=5, pady=5)
bottomEntry = Tkinter.Spinbox(barrierFrame, width=4, from_=0, to=Nx)
bottomEntry.grid(column=1, row=1, sticky='nsw', padx=5, pady=5)
ttk.Label(barrierFrame, text="Top:").grid(column=0, row=2, sticky='nes', padx=5, pady=5)
topEntry = Tkinter.Spinbox(barrierFrame, width=4, from_=0, to=Nx)
topEntry.grid(column=1, row=2, sticky='nsw', padx=5, pady=5)

gaps = [[50,70],[230,250]]
barrierFrames = []

redrawBarrierFrame()

#almost done
redrawCanvases();

#set everything in motion
root.mainloop()