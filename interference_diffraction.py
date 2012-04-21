import Tkinter, tkFileDialog, ttk, tkMessageBox
import csv
import time
#import numpy as np
from numpy import *
import numpy
from PIL import Image, ImageTk, ImageDraw, ImageFont

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

nStable = int((barrierX + sqrt((Nx-barrierX)**2 + Ny**2))*d/c/dt) #time until stability is reached: to barrier then longest way across right domain

#the wave:
lamb  = 20*d		#wavelength -- 20 points per wavelength yeilds good results
k     = 2*pi/lamb	#wavenumber
omega = 2*pi*c/lamb	#angular frequency
tau   = lamb/c		#time period

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
      redrawBarrierFrame()
      conditionalRedraw()
  
def redrawBarrierFrame():
  global barrierFrames
  global strVars
  global distances
  
  gaps.sort()
  for b in barrierFrames:
    b.destroy() #get rid of old frames
  barrierFrames = []
  strVars       = [] #need to save StringVars or else they get garbage collected
  distances     = []
  r = 3 #rows 0,1,2 already taken by widgets for adding an opeing
    
  for gap in gaps:
    bn = r-3 #barrier number
    frame = ttk.Labelframe(barrierFrame, text="Opening {}".format(bn))
    frame.grid(column=0, row=r, columnspan=2, sticky='nesw', padx=5, pady=5)
      
    top = Tkinter.IntVar()
    bottom = Tkinter.IntVar()
    distance = Tkinter.StringVar()
    setDistance(gap, distance)
    ttk.Label(frame, text="Top:").grid(column=0, row=0, sticky='nes', padx=5, pady=5)
    #having the following work is kind of tricky; the default parameter in the lambda is critical. See <http://mail.python.org/pipermail/tutor/2005-November/043360.html>
    entry = Tkinter.Spinbox(frame, width=4, textvariable=top, from_=0, to=Nx, command=lambda n=bn, tv=top: updateBarrierTop(n,tv))
    entry.bind("<Return>",lambda arg, n=bn, tv=top: updateBarrierTop(n,tv))
    entry.grid(column=1, row=0, sticky='nsw', padx=5, pady=5)
    ttk.Label(frame, text="Bottom:").grid(column=0, row=1, sticky='nes', padx=5, pady=5)
    entry = Tkinter.Spinbox(frame, width=4, textvariable=bottom, from_=0, to=Nx, command=lambda n=bn, tv=bottom: updateBarrierBottom(n,tv))
    entry.bind("<Return>",lambda arg, n=bn, tv=bottom: updateBarrierBottom(n,tv))
    entry.grid(column=1, row=1, sticky='nsw', padx=5, pady=5)
    ttk.Label(frame, textvariable=distance).grid(column=0, row=2, sticky='nes', padx=5, pady=5)
    ttk.Button(frame, text='Remove', command=lambda n=bn: removeBarrier(n)).grid(column=0, row=3, sticky='nsew', columnspan=2, padx=5, pady=5)
    top.set(gap[1])
    bottom.set(gap[0])
    strVars.append(top)
    strVars.append(bottom) 
    distances.append(distance)
      
    barrierFrames.append(frame)
    r = r + 1
      
def updateBarrierTop(barrierNumber, intVar):
  global gaps
  
  value = intVar.get()
  if ((barrierNumber == (len(gaps)-1)) and (value < Ny) and (value > gaps[-1][0])) or ((barrierNumber < (len(gaps)-1)) and (value < gaps[barrierNumber+1][0]) and (value > gaps[barrierNumber][0])):
    gaps[barrierNumber][1] = value
    conditionalRedraw()
  else:
    intVar.set(gaps[barrierNumber][1])

def updateBarrierBottom(barrierNumber, intVar):
  global gaps

  value = intVar.get()
  if ((barrierNumber == 0) and (value > 0) and (value < gaps[0][1])) or ((barrierNumber > 0) and (value > gaps[barrierNumber-1][1]) and (value < gaps[barrierNumber][1])):
    gaps[barrierNumber][0] = value
    conditionalRedraw()
  else:
    intVar.set(gaps[barrierNumber][0])    
  
def setDistance(gap, textVar):
  textVar.set(str(int(round(sqrt(((gap[0]+gap[1])/2.0-sliceY)**2+(barrierX-sliceX)**2)))))
  
def removeBarrier(barrierNumber):
  global gaps

  del gaps[barrierNumber]
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
  
  if avgSetting.get() == 'sq':
    data = 256*(float32(transpose(EzSQ)/nAveraging/maxRMSY**2))  
  else:
    data = 256*(float32(transpose(EzRMS)/maxRMSY))    
  im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
  if n <= nStable and running:
    draw = ImageDraw.Draw(im)
    draw.text((Nx/3, Ny/3), str(nStable - n), font=font, fill=255.0)
  ezRMSPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas

def applyAveraging(xS, yS, invert=False, othercoord=None):
  if avgSetting.get() == 'amp': #want to display amplitude = EzRMS*sqrt(2) with zero in middle of plot
    v = int_(numpy.round((1+EzRMS[xS,yS]*sqrt(2)/maxY)*50))
    if othercoord != None: #in this case, let's plot +/-amplitude
      othercoord.extend(othercoord[::-1]) #extend the coordinates so they're like [0,1...,Nx-1,Nx,Nx,Nx-1,...,1,0]
      v = concatenate((v,100-v[::-1])) #extend the values so they're like [V0,...,Vn,-Vn,...,-V0]
  elif avgSetting.get() == 'rms': #want to display EzRMS with zero at bottom of plot
    v = int_(numpy.round((EzRMS[xS,yS]/maxRMSY)*99))
  else: #want to display EzRMS^2 with zero at bottom of plot
    v = int_(numpy.round((EzSQ[xS,yS]/nAveraging/maxRMSY**2)*99))
    
  if invert: #mirror results across axis
    v = 100 - v
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
  y = applyAveraging(EZx_all, sliceY, invert=True, othercoord=x)
  plot(draw, x, y, 'green')
  #plot Ez
  y = int_(numpy.round((-Ez[:,sliceY,0]/maxY+1)*50))
  plot(draw, x, y, 'yellow')
  #turn plot in to a format the canvas can use
  horizPlot = ImageTk.PhotoImage(image=im)

def updateVertPlot():
  global vertPlot
  
  im = Image.new('RGB', (plotD, Ny))
  draw = ImageDraw.Draw(im)
  y = range(0,Ny)
  #plot EzRMS
  x = applyAveraging(sliceX, EZy_all, othercoord=y)
  plot(draw, x, y, 'green')
  #plot Ez
  x = int_(numpy.round((Ez[sliceX,:,0]/maxY+1)*50))
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
  invGaps = [ypos for pair in gaps for ypos in pair] #flatten gaps
  invGaps.insert(0,0) #put zero for the first item
  invGaps.append(Ny) #put Ny for the last item
  #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
  for i in range(0,len(invGaps),2):
    Ezcanvas.create_line([(barrierX,invGaps[i]),(barrierX,invGaps[i+1])], width=1, fill='red')    
    EzRMScanvas.create_line([(barrierX,invGaps[i]),(barrierX,invGaps[i+1])], width=1, fill='red')
      
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
    
  #finially, show the corridnates of the slices
  xStringVar.set("x=" + str(sliceX) + "d")
  yStringVar.set("y=" + str(sliceY) + "d")
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
  if not running:
    redrawCanvases()
      
def setSliceY(y):
  global sliceY
  
  if (y >= 0) and (y < Ny):
    sliceY = y
    for gap, dist in zip(gaps, distances):
      setDistance(gap, dist)
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
    for gap, dist in zip(gaps, distances):
      setDistance(gap, dist)
    conditionalRedraw()  
    
def resetIntensity():
  global nAveraging
  global EzSQ
  
  nAveraging = 0
  EzSQ = zeros((NEZx, NEZy))
  conditionalRedraw()
      
def reset():
  global Ez
  global t
  global maxY
    
  n = 0
  resetIntensity()
  Ez = zeros((NEZx, NEZy, 3))
  maxY = 1.0
  conditionalRedraw()
  
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
  
  nEnd = n + nStable
  fastForwarding = True
  root.after(1,fastForwardStep)
    
def fastForwardStep():
  global ezPlot
  global fastForwarding

  if n <= nEnd:
    Ezcanvas.delete('all')
    EzRMScanvas.delete('all')
    im = Image.new('RGB', (Nx,Ny))
    draw = ImageDraw.Draw(im)
    draw.text((Nx/3, Ny/3), str(nEnd-n), font=font, fill='red')
    ezPlot = ImageTk.PhotoImage(image=im)
    Ezcanvas.create_image(0,0,image=ezPlot,anchor=Tkinter.NW)
    EzRMScanvas.create_image(0,0,image=ezPlot,anchor=Tkinter.NW)
    step(avg=False)
    root.after(1,fastForwardStep)
  else: #stop fast forwarding
    fastForwarding = False
    resetIntensity()
    if running:
      run()
    else:
      start()
  
def run():
  if running:
    timer = time.clock()
    if n == nStable:
      resetIntensity()
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
  
  #next, take care of Hx and Hy using the standard Yee algorithm
  Hx[HXx_range, HXy_range] = Hx[HXx_range, HXy_range] + Db*(Ez[HXx_range, HXy_range, 0] - Ez[HXx_range, 1:(NHXy+1), 0]) #HXy_range+1
  Hy[HYx_range, HYy_range] = Hy[HYx_range, HYy_range] + Db*(Ez[1:(NHYx+1), HYy_range, 0] - Ez[HYx_range, HYy_range, 0]) #HYx_range+1
    
  #move old Ez back in the array
  Ez = roll(Ez,1,axis=2)
    
  #do the normal Yee updates on Ez in the relevant range
  Ez[EZx_range, EZy_range, 0] = Ez[EZx_range, EZy_range, 1] + Cb*(Hy[EZx_range, EZy_range] - Hy[barrierX:(NEZx-2),EZy_range] + Hx[EZx_range, 0:NEZy-2] - Hx[EZx_range, EZy_range])
    
  ##now take care of the Mur RBCs
  ##for x=NEZx-1
  rng   = slice(1,NEZy-1) #range of everything in y except the corners
  rngp1 = slice(2,NEZy)
  rngm1 = slice(0,NEZy-2)
  Ez[-1,rng,0]  = -Ez[-2,rng,2] + Ma*(Ez[-2,rng,0] + Ez[-1,rng,2]) + Mb*(Ez[-1,rng,1] + Ez[-2,rng,1]) + Mc*(Ez[-1,rngp1,1] - 2*Ez[-1,rng,1] + Ez[-1,rngm1,1] + Ez[-2,rngp1,1] - 2*Ez[-2,rng,1] + Ez[-2,rngm1,1])

  #for y=0
  rng = slice(barrierX+2,NEZx-1) #range of everything in x except right corner
  rngp1 = slice(barrierX+3,NEZx)
  rngm1 = slice(barrierX+1,NEZx-2)
  Ez[rng,0,0]  = -Ez[rng,1,2] + Ma*(Ez[rng,1,0] + Ez[rng,0,2]) + Mb*(Ez[rng,0,1] + Ez[rng,1,1]) + Mc*(Ez[rngp1,0,1] - 2*Ez[rng,0,1] + Ez[rngm1,0,1] + Ez[rngp1,1,1] - 2*Ez[rng,1,1] + Ez[rngm1,1,1])

  #for y=NEZy
  Ez[rng,-1,0]  = -Ez[rng,-2,2] + Ma*(Ez[rng,-2,0] + Ez[rng,-1,2]) + Mb*(Ez[rng,-1,1] + Ez[rng,-2,1]) + Mc*(Ez[rngp1,-1,1] - 2*Ez[rng,-1,1] + Ez[rngm1,-1,1] + Ez[rngp1,-2,1] - 2*Ez[rng,-2,1] + Ez[rngm1,-2,1])

  #now for the corners
  Ez[-1,0,0] = Ez[-2,1,2] #bottom right
  Ez[-1,-1,0] = Ez[-2,-2,2] #top right
    
  ##finially, update the time and the sources
  n = n + 1
  nAveraging = nAveraging + 1

  #update calculated part of Ez so that it displays correctly
  #the Ez source is the area [0,barrierX]
  #x and y for points on the barrier and to the left
  x, y = d * mgrid[EZx_ex_range, EZy_ex_range]
  tr = n*dt - x/c  
  
  #want wave to start gradually and propigate at speed of light
  #logistic growth makes it come in gradually and use of retarded time there and in step function at end enforces propigation
  Ez[EZx_ex_range, EZy_ex_range, 0] = sin(k*x-omega*n*dt)/(1+exp(-(tr-3*tau)/tau))*((tr > 0).astype(float))

  #now, enforce Ez=0 on barrier
  invGaps = [ypos for pair in gaps for ypos in pair] #flatten gaps
  invGaps.insert(0,0) #put zero for the first item
  invGaps.append(Ny) #put Ny for the last item
  #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
  for i in range(0,len(invGaps),2):
    Ez[barrierX, invGaps[i]:invGaps[i+1], 0] = 0     
    
  #strictly speaking, the following will blow up to infinity if you integrate forever (since the FT of a sinusoid is a delta function)
  #however, we're not that patient. plus, floating point limitations will prevent it (once the numbers are large enough, adding a small number to them won't change them)
  if avg:
    EzSQ = EzSQ + square(Ez[:,:,0])

#GLOBAL VARIABLES--------------------------------------------------------------
n = 0 #the current time step
nAveraging = 0 #number of time steps average has been running
running = False #the simulation isn't currently running
Ez = zeros((NEZx, NEZy,3)) #3rd dimension to keep track of past values of Ez.
EzSQ = zeros((NEZx, NEZy))
Hx = zeros((NHXx, NHXy))
Hy = zeros((NHYx, NHYy))
maxY = 1 #contains the largest Ez seen
font = ImageFont.truetype("BebasNeue.otf", 90)
fastForwarding = False

#THE GUI-----------------------------------------------------------------------
#The root window
root = Tkinter.Tk()
root.title("Leon's New Fangled Interference & Diffraction Simulator")

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
menubar.add_cascade(label="View", menu=viewmenu)

#the simulation menu
simmenu = Tkinter.Menu(menubar, tearoff=0)
simmenu.add_command(label="Run", accelerator="Ctrl+R", command=start)
root.bind("<Control-r>", lambda arg: start())
simmenu.add_command(label="Stop", accelerator="Ctrl+S", command=stop)
root.bind("<Control-s>", lambda arg: stop())
simmenu.add_command(label="Reset Simulation", accelerator="Ctrl+Shift+R", command=reset)
root.bind("<Control-R>", lambda arg: reset())
simmenu.add_command(label="Reset Average", accelerator="Ctrl+Shift+A", command=resetIntensity)
root.bind("<Shift-A>", lambda arg: resetIntensity())
simmenu.add_command(label="Fast Forward", accelerator="Ctrl+F", command=fastForward)
root.bind("<Control-f>", lambda arg: fastForward())
menubar.add_cascade(label="Simulation", menu=simmenu)

root.config(menu=menubar)

#The view frame
viewFrame = ttk.Labelframe(root, text='View')
viewFrame.grid(column=0, row=0, sticky='nsew',padx=5,pady=5)

Ezcanvas = Tkinter.Canvas(viewFrame, width=Nx, height=Ny)
Ezcanvas.grid(column=0, row=0, columnspan=3, rowspan=3, sticky='nsew', padx=5, pady=5)    

HorizPlotCanvas = Tkinter.Canvas(viewFrame, width=Nx, height=plotD)
HorizPlotCanvas.grid(column=0, row=3, columnspan=3, rowspan=5, sticky='nsew', padx=5, pady=5)

xStringVar = Tkinter.StringVar()
yStringVar = Tkinter.StringVar()
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

sliceY = 60 #position of horizontal slice
sliceX = Nx/2

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

#The barrier frame frame and intial barrier setup
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